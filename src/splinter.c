/**
 * SPDX-License-Identifier: GPL-2.0+
 * @file splinter.c
 * @brief Spline interpolation
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

#include "splinter.h"
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

// ********************** boundary condition **********************************

/// \brief Boundary handling function for constant extension
/// \param N the data length
/// \param i an index into the data
/// \return an index between 0 and N-1
inline static int constExt(int N, int i) {
    if(i < 0)
        return 0;
    if(i >= N)
        return N - 1;
    return i;
}

/// \brief Boundary handling function for half-sample symmetric extension
/// \param N the data length
/// \param i an index into the data
/// \return an index between 0 and N-1
inline static int hSymExt(int N, int i) {
    while(1) {
        if(i < 0)
            i = -1-i;
        else if(i >= N)
            i = (2*N-1)-i;
        else
            return i;
    }
}

/// \brief Boundary handling function for whole-sample symmetric extension
/// \param N the data length
/// \param i an index into the data
/// \return an index between 0 and N-1
inline static int wSymExt(int N, int i) {
    while(1) {
        if(i < 0)
            i = -i;
        else if(i >= N)
            i = (2*N-2)-i;
        else
            return i;
    }
}

/// \brief Boundary handling function for whole-sample symmetric extension
/// \param N the data length
/// \param i an index into the data
/// \return an index between 0 and N-1
inline static int periodicExt(int N, int i) {
    while(1) {
        if(i < 0)
            i = N+i;
        else if(i >= N)
            i = i-N;
        else
            return i;
    }
}

/// \brief Array of boundary extension methods
static int (*ExtensionMethod[4])(int, int) =
    {constExt, hSymExt, wSymExt, periodicExt};

// ********************** prefiltering exact domain ***************************

/// \brief 1D in-place exponential filter with a recursive filter pair
/// \details This is Algorithm 3 in the IPOL article.
/// \param data pointer to data to be filtered
/// \param step stride between successive elements of \a data
/// \param n number of samples of \a data
/// \param boundary the kind of boundary handling to use
/// \param alpha filter coefficient
/// \param n0 truncation index for initial values
///
/// Applies the causal recursive filter
///     1/(1 - alpha z^-1)
/// followed by the anti-causal recursive filter
///     -alpha/(1 - alpha z).
/// The coefficient alpha must satisify |alpha| < 1 for stability.
///
/// With respect to boundary handling, filtering is computed with relative
/// accuracy eps for half-and whole-sample symmetric boundaries and it
/// is exact for constant extension.  Note, however, that for constant extension
/// the infinite grid result is not exactly constant beyond the boundaries
/// (rather it decays to constant).
static void expFilter(double *data, int step, int n,
                      BoundaryExt boundary, double alpha, int n0) {
    double powAlpha=1, last=data[0];

    // avoid too large initialization
    if(n0 > n)
        n0 = n;
    if(n0 == n && boundary == BOUNDARY_WSYMMETRIC)
        n0 = n-1;
    int i, iEnd=n0*step;
    // Causal init
    switch(boundary) {
    case BOUNDARY_CONSTANT:
        last /= 1-alpha;
        break;
    case BOUNDARY_HSYMMETRIC:
        for(i=0; i<iEnd; i+=step) {
            powAlpha *= alpha;
            last += data[i]*powAlpha;
        }
        break;
    case BOUNDARY_WSYMMETRIC:
        for(i=step; i<=iEnd; i+=step) {
            powAlpha *= alpha;
            last += data[i]*powAlpha;
        }
        break;
    case BOUNDARY_PERIODIC:
        for(i=step; i<=iEnd; i+=step) {
            powAlpha *= alpha;
            last += data[step*n-i]*powAlpha;
        }
        break;
    default: assert(0); // Should never go here
        break;
    }
    data[0] = last;

    // Causal filter
    iEnd = (n-1)*step;
    for(i=step; i<iEnd; i+=step) {
        data[i] += alpha*last;
        last = data[i];
    }

    // Anti-causal init
    switch(boundary) {
    case BOUNDARY_CONSTANT:
        data[iEnd] = last = (alpha*(-data[iEnd] + (alpha - 1)*alpha*last))
            /((alpha - 1)*(alpha*alpha - 1));
        break;
    case BOUNDARY_HSYMMETRIC:
        data[iEnd] += alpha*last;
        last = data[iEnd] *= alpha/(alpha - 1);
        break;
    case BOUNDARY_WSYMMETRIC:
        data[iEnd] += alpha*last;
        data[iEnd] = last = (alpha/(alpha*alpha - 1))
            * ( data[iEnd] + alpha*data[iEnd - step] );
        break;
    case BOUNDARY_PERIODIC:
        data[iEnd] += alpha*last;
        last = data[iEnd];
        powAlpha = 1;
        for(i=0; i<n0*step; i+=step) {
            powAlpha *= alpha;
            last += data[i]*powAlpha;
        }
        data[iEnd] = last *= -alpha;
        break;
    }
    // Anti-causal filter
    for(i=iEnd-step; i>=0; i-=step) {
        data[i] = alpha*(last - data[i]);
        last = data[i];
    }
}

/// \brief Apply a cascade of exponential filters to an image
/// \details This is Algorithm 5 in the IPOL article.
/// \param data the image data
/// \param w,h image dimensions
/// \param boundary the kind of boundary handling to use
/// \param m structure with poles and number of poles
/// \param truncation array of truncation values in the initializations
static void prefiltering(double* data, int w, int h, BoundaryExt boundary,
                         const prefilter_t* m, const int* truncation) {
    int x, y, k;

    // Prefiltering of the columns
    for(x = 0; x < w; x++)
        for(k = 0; k < m->nPoles; k++)
            expFilter(data + x, w, h, boundary, m->poles[k], truncation[k]);

    // Prefiltering of the rows
    for(y = 0; y < h; y++)
        for(k = 0; k < m->nPoles; k++)
            expFilter(data+w*y, 1, w, boundary, m->poles[k], truncation[k]);

    // Normalization, twice because 2D
    if(m->normalization != 1) {
        unsigned long long factor = m->normalization*m->normalization;
        for(k = 0; k < w*h; k++)
            data[k] *= factor;
    }
}

/// \brief 1D in-place exp filter with a recursive filter pair (larger domain)
/// \details This is Algorithm 3 in the IPOL article.
/// \param data pointer to data to be filtered
/// \param step stride between successive elements of \a data
/// \param n number of samples of \a data
/// \param alpha filter coefficient
/// \param n0 truncation index for initial values
static void expFilterExt(double *data, int step, int n, double alpha, int n0) {
    int i, iIni = n0*step, iEnd = (n-1-n0)*step;

    // Initialisation at point n0 using the n0 first values
    double powAlpha=1, last=data[iIni];
    for(i=iIni-step; i>=0; i-=step) {
      powAlpha *= alpha;
      last += powAlpha*data[i];
    }
    data[iIni] = last;

    // Computation for the anti-causal initialization at n-1-n0
    double sum = 0;
    powAlpha = 1;
    for(i=iEnd+step; i<=(n-1)*step; i+=step) {
      powAlpha *= alpha;
      sum += powAlpha*data[i];
    }

    // Causal filtering from n0 to iEnd
    for(i=iIni+step; i<=iEnd; i+=step) {
        data[i] += alpha*last;
        last = data[i];
    }

    // Initialization at point n-1-n0
    last = data[iEnd] = alpha/(alpha*alpha-1)*(last+sum);

    // Anti-causal filtering
    for(i=iEnd-step; i>=iIni; i-=step) {
        data[i] = alpha*(last - data[i]);
        last = data[i];
    }
}

/// \brief Apply a cascade of exponential filters to an image (larger domain)
/// \details This is Algorithm 4 in the IPOL article.
/// \param data the image data
/// \param w,h image dimensions
/// \param boundary the kind of boundary handling to use
/// \param m structure with poles and number of poles
/// \param truncation array of truncation values in the initializations
/// \param Lprecision array of larger domain extensions
static void prefilteringExt(double* prefilt, const double* data,int w,int h,
                            BoundaryExt boundary,
                            const prefilter_t* m, const int* truncation,
                            const int* Lprecision) {
    int k, x, y;
    int nPoles = m->nPoles;
    // extended domain sizes
    int L2 = Lprecision[0];
    int w2 = w+2*L2;
    int h2 = h+2*L2;

    // extend the input data
    int (*Extension)(int, int) = ExtensionMethod[boundary];
    for(y=0; y<h2; y++) {
        int y0 = ((0<=y-L2 && y-L2<h)? y-L2: Extension(h,y-L2));
        int offset0 = w*y0;
        int offset = w2*y;
        for(x=0; x<w2; x++){
            int x0 = ((0<=x-L2 && x-L2<w)? x-L2: Extension(w,x-L2));
            prefilt[x+offset] = data[x0+offset0];
      }
    }

    // L2-Lprecision[k] = sum_{i=0}^{k-1} truncation[i] is the length of values
    // that are not used for computing the k-th application of exp filter
    if(nPoles > 0) { // security check
        // prefiltering of the columns
        for(x = 0; x < w2; x++)
            for(k = 0; k < nPoles; k++)
                expFilterExt(prefilt+x+(L2-Lprecision[k])*w2, w2,
                             h2-2*(L2-Lprecision[k]),
                             m->poles[k], truncation[k]);

        // prefiltering of the rows, needs to be computed only from
        // L3 = sum(truncation[i]) to h2-L3
        int L3 = L2-Lprecision[nPoles];
        for(y=L3; y < h2-L3; y++)
            for(k = 0; k < nPoles; k++)
                expFilterExt(prefilt + w2*y + (L2-Lprecision[k]), 1,
                             w2-2*(L2-Lprecision[k]),
                             m->poles[k], truncation[k]);

        // renormalization
        if(m->normalization != 1) {
            unsigned long long factor = m->normalization*m->normalization;
            for(y=L3; y < h2-L3; y++)
                for(x=L3; x < w2-L3; x++)
                    prefilt[x+w2*y] *= factor;
        }
    }
}

/// \brief Create a plan for spline interpolation.
/// \details This performs the prefiltering of the image and stores the result.
/// After usage by calls to function \ref splinter, the plan must be disposed of
/// with \ref splinter_destroy_plan.
/// \param in the input image (if color, in planar form: RR...RBB...BGG...G).
/// \param w number of pixels horizontally.
/// \param h number of pixels vertically.
/// \param c number of channels (usually 1 or 3).
/// \param order spline order
/// \param e rule of image extension.
/// \param eps precision required.
/// \param larger whether to compute in the original domain or in a larger one.
///
/// The usage pattern is:
/// \code
/// #include "splinter.h"
/// splinter_plan_t plan = splinter_plan(...); // Prefiltering
/// double pixOut[3];                          // Output values (#channels)
/// splinter(pixOut, 1.3, 2.4, plan);         // Interpolate at coords (1.3,2.4)
/// splinter_destroy_plan(plan);               // Free reserved memory
/// \endcode

splinter_plan_t splinter_plan(const double* in, int w, int h, int c,
                              int order, BoundaryExt e, double eps, int larger){
    splinter_plan_t plan = {.w=w, .h=h, .c=c, .shift=0};
    prefilter_t prefilter;
    plan.bspline = malloc(sizeof(Bspline));
    get_bspline(order, &prefilter, plan.bspline);

    // compute the truncation values
    int tn = prefilter.nPoles;
    int* truncation = malloc(tn*sizeof*truncation);
    if(tn > 0)
        compute_truncation(truncation, prefilter.poles, tn, eps);
    int* Lprecision=NULL;

    if(larger) {
        Lprecision = malloc((tn+1)*sizeof*Lprecision);
        Lprecision[tn] = tn;
        for(int i=tn-1; i>=0; i--)
            Lprecision[i] = Lprecision[i+1] + truncation[i];
        plan.shift = Lprecision[0];
        plan.w += 2*plan.shift;
        plan.h += 2*plan.shift;
    }

    plan.prefilt = malloc(plan.w*plan.h*c*sizeof*plan.prefilt);
    if(! larger)
        memcpy(plan.prefilt, in, w*h*c*sizeof(double));
    for(int l=0; l<c; l++) {
        if(larger)
            prefilteringExt(plan.prefilt+l*plan.w*plan.h, in+l*w*h, w, h,
                            e, &prefilter, truncation, Lprecision);
        else
            prefiltering(plan.prefilt+l*plan.w*plan.h, w, h,
                         e, &prefilter, truncation);
    }
    if(order > MAX_TABULATED_ORDER)
        free(prefilter.poles);

    plan.ext = ExtensionMethod[e];
    // B-spline of order 0 does not vanish at its support bounds
    int kWidth = (order==0)? 2: order+1;
    plan.xBuf = malloc(kWidth*sizeof*plan.xBuf);
    plan.yBuf = malloc(kWidth*sizeof*plan.yBuf);

    free(Lprecision);
    free(truncation);
    return plan;
}

/// \brief Dispose of a plan created with \ref splinter_plan.
/// \details Must be called when a plan is not used anymore.
void splinter_destroy_plan(splinter_plan_t plan) {
    if(plan.bspline->order > MAX_TABULATED_ORDER)
        free(plan.bspline->C);
    free(plan.bspline);
    free(plan.prefilt);
    free(plan.xBuf);
    free(plan.yBuf);
}

/// \brief Perform spline interpolation at coordinates (x,y).
/// \details The plan has to be created with \ref splinter_plan, which performs
/// prefiltering. The resulting pixel value is stored in \c out, which must
/// an array large enough to accomodate the number of channels of the image.
/// \param out the array (or pointer if single channel) where output values
/// are stored.
/// \remark Pixels outside the image receive the value 0. However, uncommenting
/// a single line in the function allows extrapolation with the extension
/// specified at creation of the plan.
/// \param x,y coordinates of pixel.
/// \param plan the plan create with \ref splinter_plan.
/// \details This is Algorithm 7 in the IPOL article.
void splinter(double* out, double x, double y, splinter_plan_t plan) {
    double (*betan)(double, const Bspline*) = plan.bspline->eval;
    double radius = plan.bspline->radius;

    // B-spline of order 0 does not vanish at its support bounds
    const int kWidth = (plan.bspline->order==0)? 2: plan.bspline->order+1;

    const int shift = plan.shift;
    // Shift for handling boundary condition properly in case of extrapolation
    int shift2 = (shift-plan.bspline->tn>0)? shift-plan.bspline->tn: 0;
    x += shift;
    y += shift;

    for(int c=0; c<plan.c; c++)
        out[c]=0;
    int inside = shift<=x && x<=plan.w-1-shift && shift<=y && y<=plan.h-1-shift;
    //inside=1; // Uncomment to extrapolate
    if(! inside)
        return;
    // Evaluate the kernel
    int x0 = ceil(x-radius), y0 = ceil(y-radius);
    for(int k = 0; k < kWidth; k++)
        plan.xBuf[k] = betan(x-(x0+k), plan.bspline);
    for(int k = 0; k < kWidth; k++)
        plan.yBuf[k] = betan(y-(y0+k), plan.bspline);

    // Compute the interpolated value at (x,y)
    for(int l=0; l<kWidth; l++) {
        int iY = (shift2<=y0+l && y0+l<plan.h-shift2)?
            y0+l: plan.ext(plan.h-2*shift, y0+l-shift)+shift;
        int rowOffset = plan.w*iY;

        for(int c=0; c<plan.c; c++) {
            double s=0;
            for(int k=0; k<kWidth; k++) {
                int iX = (shift2<=x0+k && x0+k<plan.w-shift2)?
                    x0+k: plan.ext(plan.w-2*shift, x0+k-shift)+shift;
                s += plan.prefilt[iX+rowOffset]*plan.xBuf[k];
            }
            out[c] += s*plan.yBuf[l];
            rowOffset += plan.w*plan.h;
        }
    }
}
