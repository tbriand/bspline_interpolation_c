/**
 * SPDX-License-Identifier: GPL-2.0+
 * @file bspline.c
 * @brief Bspline function comptutation
 * @author Thibaud Briand <thibaud.briand@enpc.fr>
 *         Pascal Monasse <monasse@imagine.enpc.fr>
 *
 * Copyright (c) 2017-2018, Thibaud Briand, Pascal Monasse
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * You should have received a copy of the GNU General Pulic License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "bspline.h"

#ifdef GSL_SUPPORT
#include <gsl/gsl_poly.h>
#endif

// ********************** coefficients and poles computation ******************

/// \brief Coefficients of the z-transform of the B-spline of order n.
/// \details These are built from the bn, values of the B-spline at integers:
/// \f[ b_i = \beta^{(n)}(i). \f]
/// This is Algorithm 9 in the IPOL article.
/// \param[out] coeff array of size n+1 (n even) or n (n odd) coefficients.
/// \param n the order of the spline.
void compute_ztrans_coeff(double* coeff, int n) {
    int tn = n/2; // tilde n

    // computation of the polynomial coefficients bn
    double *bn = malloc((tn+2)*sizeof*bn);
    double *dn = malloc((tn+1)*sizeof*dn);

    bn[0] = 1;
    dn[0] = 0.5;
    for(int k=1; k<=tn; k++) {
        bn[k] = 0;
        dn[k] = 0;
    }
    bn[tn+1] = 0;

    if(n>1) {
        double *bn_tmp = malloc((tn+1)*sizeof*bn_tmp);
        double *dn_tmp = malloc((tn+1)*sizeof*dn_tmp);
        for(int m=1; m<=n-1; m++) {
            int tm = m/2;
            double invm = 1.0/m;
            for(int k=0; k<=tm; k++) {
                bn_tmp[k] = invm*(((m+1)*0.5+k)*dn[k]+
                                  ((m+1)*0.5-k)*dn[(k>0)? k-1 :0]);
                dn_tmp[k] = invm*(((m+2)*0.5+k)*bn[k+1]+
                                   (m   *0.5-k)*bn[k]);
            }
            for(int k=0; k<=tm; k++) {
                bn[k] = bn_tmp[k];
                dn[k] = dn_tmp[k];
            }
        }
        free(bn_tmp);
        free(dn_tmp);
    }

    if(n>0) {
        double invn = 1.0/n;
        for(int k=0; k<=tn; k++)
            bn[k] = invn*(((n+1)*0.5+k)*dn[k]+
                          ((n+1)*0.5-k)*dn[(k>0)? k-1 :0]);
    }

    for(int k=0; k<=tn; k++)
        coeff[k] = bn[tn - k];
    for(int k=tn+1; k<=2*tn; k++)
        coeff[k] = coeff[2*tn-k];

    free(bn);
    free(dn);
}

#ifdef GSL_SUPPORT // Need GSL to find roots of polynomials

/// \brief Computation of the poles of B-spline interpolation.
/// \details Computation of the roots of \f$B^{(n)}\f$ or equivalently:
/// \f[z^{\tilde n} B^{(n)}(z) = b^n_{\tilde n} \left(
/// 1 + b^n_{\tilde n-1}/b^n_{\tilde n} z +\dots+
/// b^n_{0}/b^n_{\tilde n} z^{\tilde n} +\dots+
/// b^n_{\tilde n-1}/b^n_{\tilde n} z^{2 \tilde n-1} + z^{2\tilde n}\right)\f]
/// This function requires GSL support.
void compute_poles(double* poles, const double* coeff, int n) {
    int tn = n/2; // number of poles
    int count = 0; // count of roots

    // find roots
    if(n>1) {
        gsl_poly_complex_workspace * w=gsl_poly_complex_workspace_alloc(2*tn+1);
        double* z = malloc(2*(2*tn)*sizeof*z);
        gsl_poly_complex_solve(coeff, 2*tn+1, w, z);
        for(int i=2*tn-1; i>=0; i--)
            if(z[2*i]>-1.0) {
                poles[count] = z[2*i];
                count++;
            }
        gsl_poly_complex_workspace_free (w);
        free(z);
    }
    assert(count==tn);
}

#else

/// \brief Computation of the poles of B-spline interpolation.
/// \details Computation of the roots of \f$B^{(n)}\f$ or equivalently:
/// \f[z^{\tilde n} B^{(n)}(z) = b^n_{\tilde n} \left(
/// 1 + b^n_{\tilde n-1}/b^n_{\tilde n} z +\dots+
/// b^n_{0}/b^n_{\tilde n} z^{\tilde n} +\dots+
/// b^n_{\tilde n-1}/b^n_{\tilde n} z^{2 \tilde n-1} + z^{2\tilde n}\right)\f]
/// This function requires GSL support.
void compute_poles(double* poles, const double* coeff, int n) {
    (void)poles; (void)coeff; (void)n;
    exit(1);
}

#endif

/// \brief Fill \a C with polynomial coefficients of Bspline of order \a n.
/// \details This is Algorithm 2 of the IPOL article.
///
/// There are tn+1=int(n/2)+1 intervals, and a polynomial expression of degree n
/// in each, so n+1 coefficients: C is a (tn+1)*(n+1) matrix, row-major storage.
/// Interval number i is [(n+1)/2-i-1,(n+1)/2-i] (i starting at 0).
/// C_ij is the coefficient of y^j for x in the interval [(n+1)/2-i-1,(n+1)/2-i]
/// (i and j starting at 0), with y=(n+1)/2-i-x.
/// The exception is the last row C_tn,j: odd degrees do not appear (except n,
/// if n is odd), so only tn+1 (even n) or tn+2 (odd n) packed coefficients.
void compute_bspline_poly(double* C, int n) {
    int i, j, k, sign;
    long long sum;
    const int tn=n/2;

    // (1) Auxiliary tables

    // nCj[j] = n choose j; j=0:n
    long long *nCj = malloc((n+1)*sizeof*nCj);
    nCj[0]=1;
    for(j=1; j<=tn; j++)
        nCj[j] = (n-j+1)*nCj[j-1]/j;
    for(j=tn+1; j<=n; j++) // Completion by symmetry
        nCj[j] = nCj[n-j];

    // n1Ci[i] = n+1 choose i; i=0:tn
    long long *n1Ci = malloc((tn+1)*sizeof*n1Ci);
    n1Ci[0] = 1;
    for(i=1; i<=tn; i++)
        n1Ci[i] = nCj[i-1]*(n+1)/i;

    // jPi_ij = j^i, i=0:n j=1:n+1
    int n1=n+1;
    long long *jPi = malloc(n1*n1*sizeof*jPi);
    for(j=1; j<=n1; j++) // row 0
        jPi[j-1] = 1;
    for(i=1; i<=n; i++)
        for(int j=1; j<=n1; j++)
            jPi[(j-1)+n1*i] = j*jPi[(j-1)+n1*(i-1)];

    // (2) Computation of C_k,j
    for(k=0; k<tn; k++) {
        // j-th coefficient, j<n
        for(j=0; j<n; j++) {
            sum = 0;
            sign = 1;
            for(i=0;i<=k-1;i++) {
                sum += sign*n1Ci[i]*jPi[(k-i-1)+n1*(n-j)];
                sign = -sign;
            }
            C[j+n1*k] = (double)nCj[j]*sum;
        }

        // n-th coefficient
        sum = 0;
        sign = 1;
        for(i=0; i<=k; i++) {
            sum += sign*n1Ci[i];
            sign = -sign;
        }
        C[n+n1*k] = (double)sum;
    }

    // (3) Compute D_tn,j (x near 0, k=tn)

    // j-th coefficient, all j except last one
    int nodd = n%2; // Is n odd?
    double power2 = (nodd? 0.5: 0.25);
    for(j=n-2+nodd; j>=0; j-=2) {
        sum = 0;
        sign = 1;
        for(i=0; i<=tn; i++) {
            sum += sign*n1Ci[i]*jPi[(n-2*i)+n1*(n-j)];
            sign = -sign;
        }
        C[j/2+n1*tn] = power2*nCj[j]*sum;
        power2 *= 0.25;
    }

    // Last coefficient
    sum = 0;
    sign = (nodd? -1: 1);
    for(i=0; i<=tn; i++) {
        sum += sign*n1Ci[i];
        sign = -sign;
    }
    C[(tn+nodd)+n1*tn] = (double)sum;

    free(nCj);
    free(n1Ci);
    free(jPi);
}

/// \brief Evaluate Bspline at point \a x.
/// \details Find the length-1 interval x is in, wich yields the polynomial
/// coefficients to use. Evaluate with Horner's method.
static double bsplineEval(double x, const Bspline* c) {
    x = fabs(x); // symmetry
    int j, k=ceil(c->radius-x)-1;
    double *rowC, val=0; // case |x| >= (n+1)/2

    if(k==c->tn) { // x near 0: only even-degrees (except highest) monomials
        rowC = c->C+(c->order+1)*k; // Find polynomial
        val = rowC[k];
        if(c->order%2 == 1)
            val += rowC[k+1]*x;
        x=x*x;
        for(j=k-1; j>=0 ; j--)
            val = rowC[j]+val*x;
    }
    else if(k>=0) { // case 0< (n+1)/2-k-1 <= x < (n+1)/2-k
        rowC = c->C+(c->order+1)*k; // Find polynomial
        x = c->radius - x - k;
        val = rowC[c->order];
        for(j=c->order-1; j>=0; j--)
            val = rowC[j]+val*x;
    }

    return val;
}

// ******************** Truncation indices ************************************
/// \brief Compute numbers mu_i (1<=i<=tn) involved in truncation indices.
void compute_mu(double* mu, const double* poles, int tn) {
    mu[0] = 0;
    double sum = 0;
    double log_pole = log(-poles[0]);

    for(int i=1; i<tn; i++) {
        sum += 1.0/log_pole;
        log_pole = log(-poles[i]);
        mu[i] = log_pole * sum / (1 + log_pole * sum);
    }
}

/// \brief Compute truncation indices at a given precision level.
/// \details Fills the array \a trunc of length \a tn.
/// \param[out] trunc array of truncation indices
/// \param poles poles of the spline
/// \param tn tilde n (\f$\tilde n\f$), number of poles
/// \param eps desired precision level
void compute_truncation(int* trunc, const double* poles, int tn, double eps) {
    // weight computations
    double *mu = malloc(tn * sizeof * mu);
    compute_mu(mu, poles, tn);

    // compute \rho^{(n)}
    double rhon = 1.0;
    for(int i=0; i<tn; i++)
        rhon *= (1 + poles[i]) / (1 - poles[i]);
    rhon *= rhon;

    // compute truncation indices
    double log_eps_cte = log(eps * rhon * rhon * 0.5);
    double alpha;
    double prod_mu = 1;
    for(int i=tn-1; i>=0; i--) {
        alpha = poles[i];
        trunc[i] = 1+ floor( (log_eps_cte + log((1-alpha)*(1-mu[i])*prod_mu))
                             /log(-alpha) );
        prod_mu *= mu[i];
    }

    free(mu);
}

// ********************* Tabulation of splines of small orders ****************

/// Constant B-spline function (KernelRadius = 0.5)
static double BSpline0(double x, const Bspline *c) {
    (void)c; // Shut up compiler warning
    x = fabs(x);
    if(x < 0.5)
        return 1;
    if(x==0.5)
        return 0.5;
    return 0;
}

/// Linear B-spline function (KernelRadius = 1)
static double BSpline1(double x, const Bspline *c) {
    (void)c; // Shut up compiler warning
    x = fabs(x);
    if(x < 1)
        return 1-x;
    return 0;
}

/// Quadratic B-spline function (KernelRadius = 1.5)
static double BSpline2Poles[1] =
    {-1.715728752538099e-1}; // -3+sqrt(8)
static double BSpline2(double x, const Bspline *c) {
    (void)c; // Shut up compiler warning
    x = fabs(x);
    if(x < 0.5)
        return 1.5 - 2*x*x;
    if(x < 1.5) {
        x = 1.5 - x;
        return x*x;
    }
    return 0;
}

/// Cubic B-spline function (KernelRadius = 2)
static double BSpline3Poles[1] =
    {-2.679491924311227e-1}; // -2+sqrt(3)
static double BSpline3(double x, const Bspline *c) {
    (void)c; // Shut up compiler warning
    x = fabs(x);
    if(x < 1)
        return (4 + (-6 + 3*x)*x*x);
    if(x < 2) {
        x = 2 - x;
        return x*x*x;
    }
    return 0;
}

/// Quartic B-spline function (KernelRadius = 2.5)
static double BSpline4Poles[2] =
    {-3.6134122590021989e-1,-1.3725429297339109e-2};
static double BSpline4(double x, const Bspline *c) {
    (void)c; // Shut up compiler warning
    x = fabs(x);
    if(x <= 0.5) {
        x *= x;
        return (14.375 + (-15 + 6*x)*x);
    }
    if(x < 1.5) {
        x = 1.5 - x;
        return (1 + (4 + (6 + (4 - 4*x)*x)*x)*x);
    }
    if(x < 2.5) {
        x = 2.5 - x;
        x *= x;
        return x*x;
    }
    return 0;
}

/// Quintic B-spline function (KernelRadius = 3)
static double BSpline5Poles[2] =
    {-4.305753470999738e-1,  // sqrt(105)/2+sqrt(135-13*sqrt(105))/sqrt(2)-13/2.
     -4.309628820326465e-2}; // sqrt(13*sqrt(105)+135)/sqrt(2)-sqrt(105)/2-13/2.
static double BSpline5(double x, const Bspline *c) {
    (void)c; // Shut up compiler warning
    x = fabs(x);
    if(x <= 1) {
        double x2 = x*x;
        return (((-10*x + 30)*x2 - 60)*x2 + 66);
    }
    if(x < 2) {
        x = 2 - x;
        return (1 + (5 + (10 + (10 + (5 - 5*x)*x)*x)*x)*x);
    }
    if(x < 3) {
        x = 3 - x;
        double x2 = x*x;
        return x2*x2*x;
    }
    return 0;
}

/// Sextic B-spline function (KernelRadius = 3.5)
static double BSpline6Poles[3] =
    {-4.8829458930303893e-1,-8.1679271076238694e-2,-1.4141518083257976e-3};
static double BSpline6(double x, const Bspline *c) {
    (void)c; // Shut up compiler warning
    x = fabs(x);
    if(x <= 0.5) {
        x *= x;
        return (367.9375 + (-288.75 + (105 - 20*x)*x)*x);
    }
    if(x < 1.5) {
        x = 1.5 - x;
        return (57 + (150 + (135 + (20 + (-45 + (-30 + 15*x)*x)*x)*x)*x)*x);
    }
    if(x < 2.5) {
        x = 2.5 - x;
        return (1 + (6 + (15 + ( 20 + (15 + (6 - 6*x)*x)*x)*x)*x)*x);
    }
    if(x < 3.5) {
        x = 3.5 - x;
        x = x*x;
        return x*x*x;
    }
    return 0;
}

/// Septic B-spline function (KernelRadius = 4)
static double BSpline7Poles[3] =
    {-5.352804307964382e-1, -1.225546151923267e-1,-9.148694809608277e-3};
static double BSpline7(double x, const Bspline *c) {
    (void)c; // Shut up compiler warning
    x = fabs(x);
    if(x <= 1) {
        double x2 = x*x;
        return ((((35*x - 140)*x2 + 560)*x2 - 1680)*x2 + 2416);
    }
    if(x < 2) {
        x = 2 - x;
        return (120 + (392 + (504 + (280 + (-84 + (-42 +
            21*x)*x)*x*x)*x)*x)*x);
    }
    if(x < 3) {
        x = 3 - x;
        return (((((((-7*x + 7)*x + 21)*x + 35)*x + 35)*x
            + 21)*x + 7)*x + 1);
    }
    if(x < 4) {
        x = 4 - x;
        double x2 = x*x;
        return x2*x2*x2*x;
    }
    return 0;
}

/// Octic B-spline function (KernelRadius = 4.5)
static double BSpline8Poles[4] =
    {-5.7468690924876376e-1,-1.6303526929728299e-1,
     -2.3632294694844336e-2,-1.5382131064168442e-4};
static double BSpline8(double x, const Bspline *c) {
    (void)c; // Shut up compiler warning
    x = fabs(x);
    if(x < 0.5) {
        x *= x;
        return (18261.7734375 + (-11379.375 + (3386.25 + (-630 + 70*x)*x)*x)*x);
    }
    if(x <= 1.5) {
        x = 1.5 - x;
        return (4293 + (8568 + (5292 + (-504 + (-1890 + (-504 + (252 + (168
                - 56*x)*x)*x)*x)*x)*x)*x)*x);
    }
    if(x <= 2.5) {
        x = 2.5 - x;
        return (247 + (952 + (1540 + (1288 + (490 + (-56 +(-140 + (-56
                + 28*x)*x)*x)*x)*x)*x)*x)*x);
    }
    if(x < 3.5) {
        x = 3.5 - x;
        return (1+ (8+ (28+ (56+ (70+ (56+ (28+ (8 - 8*x)*x)*x)*x)*x)*x)*x)*x);
    }
    if(x < 4.5) {
        x = 4.5 - x;
        x = x*x; x = x*x; // x^4
        return x*x;
    }
    return 0;
}

/// Nonic B-spline function (KernelRadius = 5)
static double BSpline9Poles[4] =
    {-6.079973891686259e-1,-2.017505201931532e-1,
     -4.322260854048175e-2,-2.121306903180818e-3};
static double BSpline9(double x, const Bspline *c) {
    (void)c; // Shut up compiler warning
    x = fabs(x);
    if(x <= 1) {
        double x2 = x*x;
        return (((((-63*x + 315)*x2 - 2100)*x2 + 11970)*x2
            - 44100)*x2 + 78095)*2;
    }
    if(x <= 2) {
        x = 2 - x;
        return (14608 + (36414 + (34272 + (11256 + (-4032 + (-4284 + (-672
                + (504 + (252 - 84*x)*x)*x)*x)*x)*x)*x)*x)*x);
    }
    if(x <= 3) {
        x = 3 - x;
        return (502 + (2214 + (4248 + (4536 + (2772 + (756 + (-168 + (-216
                + (-72 + 36*x)*x)*x)*x)*x)*x)*x)*x)*x);
    }
    if(x < 4) {
        x = 4 - x;
        return (1 + (9 + (36 + (84 + (126 + (126 + (84 + (36 + (9
                - 9*x)*x)*x)*x)*x)*x)*x)*x)*x);
    }
    if(x < 5) {
        x = 5 - x;
        double x3 = x*x*x;
        return x3*x3*x3;
    }
    return 0;
}

/// 10th-Degree B-spline function (KernelRadius = 5.5)
static double BSpline10Poles[5] =
    {-6.365506639694650e-1,-2.381827983775487e-1,-6.572703322831758e-2,
     -7.528194675547741e-3,-1.698276282327549e-5};
static double BSpline10(double x, const Bspline *c) {
    (void)c; // Shut up compiler warning
    x = fabs(x);
    if(x < 0.5) {
        x *= x;
        return (1491301.23828125 + (-769825.546875 + (191585.625 + (-30607.5
                + (3465 - 252*x)*x)*x)*x)*x);
    }
    if(x <= 1.5) {
        x = 1.5 - x;
        return (455192+ (736260+ (327600+ (-95760 + (-119280 + (-13608
                + (16800+ (5040+ (-1260+(-840+210*x)*x)*x)*x)*x)*x)*x)*x)*x)*x);
    }
    if(x <= 2.5) {
        x = 2.5 - x;
        return (47840 + (141060 + (171000 + (100080 + (16800 + (-13608 + (-8400
                + (-720 + (900 + (360 - 120*x)*x)*x)*x)*x)*x)*x)*x)*x)*x);
    }
    if(x <= 3.5) {
        x = 3.5 - x;
        return (1013 + (5010 + (11025 + (14040 + (11130 + (5292 + (1050
                + (-360 + (-315 + (-90 + 45*x)*x)*x)*x)*x)*x)*x)*x)*x)*x);
    }
    if(x < 4.5) {
        x = 4.5 - x;
        return (1 + (10 + (45 + (120 + (210 + (252 + (210 + (120 + (45 + (10
                -10*x)*x)*x)*x)*x)*x)*x)*x)*x)*x);
    }
    if(x < 5.5) {
        x = 5.5 - x;
        x = x*x; // x2;
        double x4 = x*x;
        return x4*x4*x;
    }
    return 0;
}

/// 11th-Degree B-spline function (KernelRadius = 6)
static double BSpline11Poles[5] =
    {-6.612660689007345e-1,-2.721803492947859e-1,-8.975959979371331e-2,
     -1.666962736623466e-2,-5.105575344465021e-4};
static double BSpline11(double x, const Bspline *c) {
    (void)c; // Shut up compiler warning
    x = fabs(x);
    if(x <= 1) {
        double x2 = x*x;
        return (15724248 + (-7475160 + (1718640 + (-255024 + (27720
            + (-2772 + 462*x)*x2)*x2)*x2)*x2)*x2);
    }
    if(x <= 2) {
        x = 2 - x;
        return (2203488 + (4480872 + (3273600 + (574200 + (-538560
            + (-299376 + (39600 + (7920 + (-2640 + (-1320
            + 330*x)*x)*x)*x)*x*x)*x)*x)*x)*x)*x);
    }
    if(x <= 3) {
        x = 3 - x;
        return (152637 + (515097 + (748275 + (586575 + (236610 + (12474
            + (-34650 + (-14850 + (-495 + (1485
            + (495-165*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)*x);
    }
    if(x <= 4) {
        x = 4 - x;
        return (2036 + (11132 + (27500 + (40260 + (38280 + (24024 + (9240
            + (1320 + (-660 + (-440 + (-110
            + 55*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)*x);
    }
    if(x < 5) {
        x = 5 - x;
        return (1 + (11 + (55 + (165 + (330 + (462 + (462 + (330 + (165
            + (55 + (11 - 11*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)*x);
    }
    if(x < 6) {
        x = 6 - x;
        double x2 = x*x;
        double x4 = x2*x2;
        return x4*x4*x2*x;
    }
    return 0;
}

static double (*BSplineTable[MAX_TABULATED_ORDER+1])(double, const Bspline*) = {
    BSpline0, BSpline1, BSpline2, BSpline3, BSpline4, BSpline5,
    BSpline6, BSpline7, BSpline8, BSpline9, BSpline10, BSpline11
};

static prefilter_t InterpMethodTable[MAX_TABULATED_ORDER+1] = {
    {0, NULL, 1},
    {0, NULL, 1},
    {1, BSpline2Poles,     4},
    {1, BSpline3Poles,     1},
    {2, BSpline4Poles,    16},
    {2, BSpline5Poles,     1},
    {3, BSpline6Poles,    64},
    {3, BSpline7Poles,     1},
    {4, BSpline8Poles,   256},
    {4, BSpline9Poles,     1},
    {5, BSpline10Poles, 1024},
    {5, BSpline11Poles,    1}
};

/// \brief Compute prefiltering parameters and spline coefficients.
/// \details Up to \c MAX_TABULATED_ORDER, no computation is necessary and
/// fields of \a p and \a s are not memory-allocated.
/// \param n spline order
/// \param[out] p prefiltering parameters
/// \param[out] s spline function coefficients
void get_bspline(int n, prefilter_t* p, Bspline* s) {
    Bspline tmp = {n, 0.5*(n+1), n/2, NULL, NULL};
    *s = tmp;
    if(n <= MAX_TABULATED_ORDER) {
        *p = InterpMethodTable[n];
        s->eval = BSplineTable[n];
    } else {
        s->order = n;
        s->eval = bsplineEval;
        s->radius = 0.5*(n+1);
        p->nPoles = s->tn = n/2;

        // computation of Bn and the normalization constant
        // in practice gamma_n is not used
        // the normalization constant is 2^n for n even and 1 for n odd
        double* zcoeff = malloc((2*p->nPoles+1)*sizeof*zcoeff);
        compute_ztrans_coeff(zcoeff, n);
        p->normalization = 1;
        if(n%2 == 0) // even order n: normalization = 2^n
            for(int i=0; i<n; i++)
                p->normalization *= 2;

        // Computation of the poles
        p->poles = malloc(p->nPoles*sizeof*p->poles);
#ifdef GSL_SUPPORT
        compute_poles(p->poles, zcoeff, n);
#else
        fprintf(stderr, "The program was built without GSL support: ");
        fprintf(stderr, "Use bspline order at most %i.\n", MAX_TABULATED_ORDER);
        exit(1);
#endif
        free(zcoeff);

        // Computation of the kernel
        s->C = malloc(((n+1)*s->tn+floor(s->radius)+1)*sizeof*s->C);
        compute_bspline_poly(s->C, n);
    }
}
