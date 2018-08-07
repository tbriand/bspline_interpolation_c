/**
 * SPDX-License-Identifier: GPL-2.0+
 * @file splinter_transform.c
 * @brief Apply homography using spline interpolation
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

#include "splinter_transform.h"
#include "homography_tools.h"
#include <stdlib.h>

/// Apply homography with spline interpolation to an image.
void splinter_homography(double *out, const double *in,
                         int w, int h, int c,
                         int n, BoundaryExt boundary, double eps,
                         int larger, const double H[9]) {
    // invert homography
    double iH[9];
    invert_homography(iH, H);

    splinter_plan_t plan = splinter_plan(in,w,h,c, n, boundary, eps, larger);

    // computation of the pixel locations
    double p[2], q[2];
    double* outp = malloc(c*sizeof*outp);
    for(int j = 0; j < h; j++) {
        p[1] = j;
        for(int i = 0; i < w; i++) {
            p[0] = i;
            apply_homography(q, p, iH);
            splinter(outp, q[0], q[1], plan);
            for(int k=0; k<c; k++)
                out[k*w*h] = outp[k];
            ++out;
        }
    }
    free(outp);

    splinter_destroy_plan(plan);
}
