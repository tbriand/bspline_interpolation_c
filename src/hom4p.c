/**
 * SPDX-License-Identifier: GPL-2.0+
 * @file hom4p.c
 * @brief Computation of the homography from 4 correspondances
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
#include <math.h>
#include "homography_tools.h"

/// Computation of the homography from 4 correspondances
int main(int c, char *v[])
{
    if (c != 17) {
            fprintf(stderr, "usage:\n\t%s x1 y1 hx1 hy1 x2 y2 hx2 hy2 x3 y3 hx3 hy3\n", *v);
            fprintf(stderr, "x4 y4 hx4 hy4\n");
            return EXIT_FAILURE;
    }

    // read parameters
    double x1 = atof(v[1]);
    double y1 = atof(v[2]);
    double hx1 = atof(v[3]);
    double hy1 = atof(v[4]);
    double x2 = atof(v[5]);
    double y2 = atof(v[6]);
    double hx2 = atof(v[7]);
    double hy2 = atof(v[8]);
    double x3 = atof(v[9]);
    double y3 = atof(v[10]);
    double hx3 = atof(v[11]);
    double hy3 = atof(v[12]);
    double x4 = atof(v[13]);
    double y4 = atof(v[14]);
    double hx4 = atof(v[15]);
    double hy4 = atof(v[16]);

    // declaration
    double initial[4][2] = {{x1,y1},{x2,y2},{x3,y3},{x4,y4}};
    double final[4][2] = {{hx1,hy1},{hx2,hy2},{hx3,hy3},{hx4,hy4}};
    double R[3][3];

    // computation of the homography
    homography_from_4corresp(
                    initial[0], initial[1], initial[2], initial[3],
                    final[0], final[1], final[2], final[3], R);

    // print the homography in stdout
    for (int i=0; i<9; i++)
        printf("%1.16lg%c", R[i/3][i%3], i==8?'\n':' ');

    return EXIT_SUCCESS;
}
