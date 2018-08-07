/**
 * SPDX-License-Identifier: GPL-2.0+
 * @file bspline_main.c
 * @brief Main program for image homography transform with bspline interpolation
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
#include "iio.h"
#include "splinter_transform.h"
#include "xmtime.h"

/// \mainpage Splinter: spline interpolation of images.
/// This program performs a homographic transformation of an image using
/// B-spline interpolation. The focus is in the B-spline interpolation, where
/// different boundary conditions can be assumed, such as periodic extension.
/// The prefiltering involves applying exponential filters, which are infinite
/// support. The exact initialization can be approximated up to an arbitrary
/// precision specified by the user.
///
/// The interpolation logic is fully contained in the Splinter library, which
/// has no dependency. This library consists of two source files and two
/// associated header files: bspline.[hc] and splinter.[hc]. The usage pattern
/// is inspired by the FFTW library:
/// \code
/// #include "splinter.h"
/// splinter_plan_t plan = splinter_plan(...); // Prefiltering
/// double pixOut[3];                          // Output values (#channels)
/// splinter(pixOut, 1.3, 2.4, plan);        // Interpolate at coords (1.3,2.4)
/// splinter_destroy_plan(plan);               // Free reserved memory
/// \endcode
/// Between creation of the plan and its destruction, any number of
/// interpolations can be performed by calling the \ref splinter function.
/// Multiple image channels are supported, in which case the channels are
/// supposed to be consecutive in memory: R...RG...GB...B for example.
///
/// This is the source code of an IPOL article, where all algorithmic aspects
/// are explained.

/// Print usage information
static void usage(const char* argv0) {
    //positions               0       1
    fprintf(stderr, "  Usage: %s \"homography\""
            " in out [order boundary eps larger]\n", argv0);
    //        2  3    4     5        6   7
    fprintf(stderr, "Homographic transformation of an image");
    fprintf(stderr, " using B-spline interpolation\n\n");
    fprintf(stderr, "homography: 9 matrix coefficients");
    fprintf(stderr, " (\"h11 h12 h13; h21 h22 h23; h31 h32 h33\")\n");
    fprintf(stderr, "in       : filename of the input image\n");
    fprintf(stderr, "out      : filename of the output image\n");
    fprintf(stderr, "  Options:\n");
    fprintf(stderr, "order    : order of interpolation (integer between ");
    fprintf(stderr, "0 and %i, default %i)\n", MAX_ORDER, MAX_TABULATED_ORDER);
    fprintf(stderr, "boundary : boundary extension");
    fprintf(stderr, " (constant");
    fprintf(stderr, ", periodic");
    fprintf(stderr, ", hsymmetric*");
    fprintf(stderr, ", wsymmetric)\n");
    fprintf(stderr, "eps      : relative precision (float, default 6)");
    fprintf(stderr, " (eps>=1 means 10^-eps)\n");
    fprintf(stderr, "larger   : compute on exact (0*) or larger domain (1)\n");
    fprintf(stderr, "  *default parameters\n");
}

/// Parse double values. They can be separated by spaces or a separator sign
static int parse_doubles(double *t, int nmax, const char *s) {
    int i = 0, w;
    char c;
    while(i<nmax && 1<=sscanf(s, "%lg%c%n", t+i, &c, &w)) {
        i += 1;
        s += w;
    }
    return i;
}

/// Read boundary extension
static BoundaryExt read_ext(const char* boundary, int* larger) {
    if(0 == strncmp(boundary, "constant", strlen(boundary))) {
        if (*larger == 0) {
            fprintf(stderr,"The constant extension is not compatible with ");
            fprintf(stderr,"computations in the exact domain.\n");
            fprintf(stderr,"\tParameter 'larger' changed to 1.\n");
            *larger = 1;
        }
        return BOUNDARY_CONSTANT;
    }
    if(0 == strncmp(boundary, "per", strlen(boundary)))
        return BOUNDARY_PERIODIC;
    if(0 == strncmp(boundary, "hsymmetric", strlen(boundary)))
        return BOUNDARY_HSYMMETRIC;
    if(0 == strncmp(boundary, "wsymmetric", strlen(boundary)))
        return BOUNDARY_WSYMMETRIC;
    fprintf(stderr,"Unknown boundary condition %s\n",boundary);
    exit(EXIT_FAILURE);
}

// eps <-> 10^(-eps) if eps>=1
static double fix_precision(double eps) {
    if(eps >= 1) {
        double tmp = 1;
        for(int i=0; i<eps; i++)
            tmp *= 0.1;
        eps = tmp;
    }
    return eps;
}

/// Apply homogaphy to an image using spline interpolation
int main(int argc, char *argv[]) {
    if(! (4<=argc && argc<=8)) {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    // Read parameters
    char *input_params = argv[1];
    char *filename_in = argv[2];
    char *filename_out = argv[3];
    int order = (argc>4? atoi(argv[4]): 11);
    char *boundary = (argc>5? argv[5]: "hsym");
    double eps = (argc>6? atof(argv[6]): 6);
    int larger = (argc>7? atoi(argv[7]): 0);

    if(order > MAX_ORDER) {
        fprintf(stderr,"The maximal order authorized is %i\n", MAX_ORDER);
        return EXIT_FAILURE;
    }

    eps = fix_precision(eps);
    BoundaryExt ext = read_ext(boundary, &larger);

    // Read transformation
    double homo[9];
    int nparams = parse_doubles(homo, 9, input_params);
    if(nparams != 9) {
        fprintf(stderr,"Wrong number of parameters in homography\n");
        return EXIT_FAILURE;
    }

    // Read input image
    int w, h, c;
    double *in = iio_read_image_double_split(filename_in, &w, &h, &c);

    int Npixels = w*h*c;
    double *out = malloc(Npixels*sizeof*out);

    unsigned long t0 = xmtime();
    splinter_homography(out, in, w,h,c, order, ext, eps, larger, homo);
    fprintf(stderr, "interpolation: %.3f s\n", (xmtime()-t0)/1000.0f);

    iio_write_image_double_split(filename_out, out, w, h, c);

    free(in);
    free(out);

    return EXIT_SUCCESS;
}
