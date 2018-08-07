/**
 * SPDX-License-Identifier: GPL-2.0+
 * @file compute_bspline.c
 * @brief Compute and display Bspline information
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "bspline.h"

/// Information about splines: prefiltering and polynomial expression
int main(int c, char *v[])
{
    // Display usage
    if(c < 2) {
        fprintf(stderr, "usage:\n\t%s order [eps]\n", *v);
        return EXIT_FAILURE;
    }

    // Read parameters
    int n = atoi(v[1]);
    double eps = (c > 2) ? atof(v[2]) : -1;

    int tn = n/2;
    double tn2 = (n+1)*0.5;

    // Computations
    double *ztransCoeff = malloc((2*tn+1)*sizeof*ztransCoeff);
    double *C = malloc(((n+1)*tn+floor(tn2)+1)*sizeof*C);
    double *poles = malloc(tn*sizeof*poles);
    double *mu = malloc(tn*sizeof*mu);
    int *truncation = malloc(tn*sizeof*truncation);

    int i, k;

    // Computation of Bn (z-transform of the prefiltering)
    compute_ztrans_coeff(ztransCoeff, n);

    // Computation of the kernel polynomial expression
    compute_bspline_poly(C, n);

    if(tn>0) {
        compute_poles(poles, ztransCoeff, n); // Poles of Bn
        compute_mu(mu, poles, tn); // Coefficients used in truncation

        // Computation of the truncation indices
        if(eps > 0) {
            if(eps >= 1) { // eps <-> 10^(-eps)
                double tmp = 1;
                for (int i=0; i<eps; i++)
                    tmp *= 0.1;
                eps = tmp;
            }
            compute_truncation(truncation, poles, tn, eps);
        }
    }

    printf("Coefficient of Bn:");
    for(k=0; k<2*tn+1; k++)
        printf(" %1.16lg",ztransCoeff[k]);
    printf("\n");

    printf("Piecewise polynomials:\n");
    k=tn;
    printf("* Case 0<=x<%g: ", tn2-k);
    printf("%1.16g ", C[0+(n+1)*k]);
    for(i=1; i<=tn; i++)
        printf("%+1.16gx^%i ", C[i+(n+1)*k], i*2);
    if(n%2)
        printf("%+1.16gx^%i", C[i+(n+1)*k], 2*tn+1);
    printf("\n");
    for(k=tn-1; k>=0; k--) {
        printf("* Case %g<=x<%g", tn2-k-1, tn2-k);
        printf(" (y=%g-x): ",tn2-k);
        printf("%.0f ",C[0 +(n+1)*k]);
        for(i=1; i<=n; i++)
            printf("%+.0fy^%i ", C[i+(n+1)*k], i);
        printf("\n");
    }

    if(tn>0) {
        printf("Poles:");
        for(i=0; i<tn; i++)
            printf(" %1.16g", poles[i]);
        printf("\n");

        printf("Mu coefficients:");
        for(i=0; i<tn; i++)
            printf(" %1.16g", mu[i]);
        printf("\n");

        if(eps > 0) {
            printf("Truncation indices for precision %1.16g:", eps);
            for(i=0; i<tn; i++)
                printf(" %i",truncation[i]);
            printf("\n");
        }
    }

    free(ztransCoeff);
    free(poles);
    free(C);
    free(mu);
    free(truncation);

    return EXIT_SUCCESS;
}
