/**
 * SPDX-License-Identifier: GPL-2.0+
 * @file splinter_transform.h
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

#ifndef SPLINTERTRANSFORM_H
#define SPLINTERTRANSFORM_H

#include "splinter.h"

void splinter_homography(double *out, const double *in, int w, int h, int c,
                         int order, BoundaryExt boundary, double eps,
                         int larger, const double homo[9]);

#endif
