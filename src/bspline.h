/**
 * SPDX-License-Identifier: GPL-2.0+
 * @file bspline.h
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

#ifndef BSPLINE_H
#define BSPLINE_H

#define MAX_TABULATED_ORDER 11 ///< Maximum order of tabulated splines
#define MAX_ORDER 16 ///< Max spline order with guaranty of truncation precision

/** Info for B-spline prefiltering */
typedef struct {
    int nPoles; ///< Number of poles = [order/2]
    double* poles; ///< Array of poles in increasing order
    unsigned long long normalization; ///< Normalization constant for prefilter
} prefilter_t;

/** Info for B-spline interpolation */
typedef struct Bspline_s Bspline;
/** B-spline structure, use shorter name \a Bspline */
struct Bspline_s {
    int order; ///< Spline order
    double radius; ///< Radius of B-spline support
    int tn; ///< Tilde n = [n/2]
    double* C; ///< Array of polynomial coefficients
    double (*eval)(double x, const Bspline* self); ///< Evaluation function
};

void compute_bspline_poly(double* C, int n);
void compute_ztrans_coeff(double* ztransCoeff, int n);
void compute_poles(double* poles, const double* ztransCoeff, int n);

void compute_mu(double* mu, const double* poles, int tn);
void compute_truncation(int* trunc, const double* poles, int tn, double eps);

void get_bspline(int n, prefilter_t* p, Bspline* s);

#endif
