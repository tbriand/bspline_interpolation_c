/**
 * SPDX-License-Identifier: GPL-2.0+
 * @file splinter.h
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

#ifndef SPLINTER_H
#define SPLINTER_H

#include "bspline.h"

/// Boundary extension method used in prefiltering
typedef enum {
    BOUNDARY_CONSTANT = 0,   ///< constant value
    BOUNDARY_HSYMMETRIC = 1, ///< half-symmetric
    BOUNDARY_WSYMMETRIC = 2, ///< whole-symmetric
    BOUNDARY_PERIODIC = 3    ///< periodic
} BoundaryExt;

/// \brief Opaque structure, intended to be used for spline interpolation.
/// \details The usage pattern is modeled after FFTW (http://www.fftw.org).
/// To interpolate, the user must first create a plan with \ref splinter_plan.
/// Calls to function \ref splinter, returning the interpolation value at points
/// (x,y), can then be performed.
/// At the end, disposal is achieved by \ref splinter_destroy_plan.
typedef struct {
    double* prefilt; ///< prefiltered image
    int w,h,c; ///< width,height,channels
    int shift; ///< shift in each channel
    Bspline* bspline; ///< Bspline kernel
    int (*ext)(int, int); ///< get pixels of extended image
    double *xBuf, *yBuf; ///< buffers for computation (internal usage)
} splinter_plan_t;

splinter_plan_t splinter_plan(const double* in, int w, int h, int c,
                              int order, BoundaryExt e, double eps, int larger);
void splinter_destroy_plan(splinter_plan_t plan);

void splinter(double* out, double x, double y, splinter_plan_t plan);

#endif
