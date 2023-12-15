/**
 * SPDX-License-Identifier: LGPL-3.0-or-later
 * @file homography_tools.h
 * @brief Routines for computing homography related stuff
 * @author Thibaud Briand <thibaud.briand@enpc.fr>
 *         Pascal Monasse <monasse@imagine.enpc.fr>
 *
 * Copyright (c) 2017-2023, Thibaud Briand, Pascal Monasse
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Pulic License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef HOMOGRAPHYTOOLS_H
#define HOMOGRAPHYTOOLS_H

void invert_homography(double iH[9], const double H[9]);
void apply_homography(double y[2], const double x[2], const double H[9]);
void homography_from_4corresp(const double *a, const double *b,
                              const double *c, const double *d,
                              const double *x, const double *y,
                              const double *z, const double *w,
                              double R[3][3]);

#endif
