/**
 * SPDX-License-Identifier: GPL-2.0+
 * @file homography_tools.c
 * @brief Routines for computing homography related stuff
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

/// Compute inverse homography
void invert_homography(double iH[9], const double H[9]) {
  double det = H[0] * (H[4]*H[8] - H[5] * H[7]);
  det -= H[1] * (H[3]*H[8] - H[5] * H[6]);
  det += H[2] * (H[3]*H[7] - H[4] * H[6]);

  double tmp = 1.0/det;

  iH[0] = tmp * (H[4] * H[8] - H[5] * H[7]);
  iH[3] = tmp * (H[5] * H[6] - H[3] * H[8]);
  iH[6] = tmp * (H[3] * H[7] - H[4] * H[6]);

  iH[1] = tmp * (H[2] * H[7] - H[1] * H[8]);
  iH[4] = tmp * (H[0] * H[8] - H[2] * H[6]);
  iH[7] = tmp * (H[1] * H[6] - H[0] * H[7]);

  iH[2] = tmp * (H[1] * H[5] - H[2] * H[4]);
  iH[5] = tmp * (H[2] * H[3] - H[0] * H[5]);
  iH[8] = tmp * (H[0] * H[4] - H[1] * H[3]);
}

/// Apply homography to a 2D point.
void apply_homography(double y[2], const double x[2], const double H[9]) {
    double z = 1.0/(H[6]*x[0] + H[7]*x[1] + H[8]);
    y[0] = z*(H[0]*x[0] + H[1]*x[1] + H[2]);
    y[1] = z*(H[3]*x[0] + H[4]*x[1] + H[5]);
}

/// Compute the homography sending [0,0] , [0,1], [1,1] and [1,0] to x,y,z,w.
static void homography_from_4pt(const double *x, const double *y,
                                const double *z, const double *w,
                                double cgret[8]) {
    double t1 = x[0];
    double t2 = z[0];
    double t4 = y[1];
    double t5 = t1 * t2 * t4;
    double t6 = w[1];
    double t7 = t1 * t6;
    double t8 = t2 * t7;
    double t9 = z[1];
    double t10 = t1 * t9;
    double t11 = y[0];
    double t14 = x[1];
    double t15 = w[0];
    double t16 = t14 * t15;
    double t18 = t16 * t11;
    double t20 = t15 * t11 * t9;
    double t21 = t15 * t4;
    double t24 = t15 * t9;
    double t25 = t2 * t4;
    double t26 = t6 * t2;
    double t27 = t6 * t11;
    double t28 = t9 * t11;
    double t30 = 0.1e1 / (-t24 + t21 - t25 + t26 - t27 + t28);
    double t32 = t1 * t15;
    double t35 = t14 * t11;
    double t41 = t4 * t1;
    double t42 = t6 * t41;
    double t43 = t14 * t2;
    double t46 = t16 * t9;
    double t48 = t14 * t9 * t11;
    double t51 = t4 * t6 * t2;
    double t55 = t6 * t14;
    cgret[0] = -(-t5 + t8 + t10 * t11 - t11 * t7 - t16 * t2 + t18 - t20 + t21 * t2) * t30;
    cgret[1] = (t5 - t8 - t32 * t4 + t32 * t9 + t18 - t2 * t35 + t27 * t2 - t20) * t30;
    cgret[2] = t1;
    cgret[3] = (-t9 * t7 + t42 + t43 * t4 - t16 * t4 + t46 - t48 + t27 * t9 - t51) * t30;
    cgret[4] = (-t42 + t41 * t9 - t55 * t2 + t46 - t48 + t55 * t11 + t51 - t21 * t9) * t30;
    cgret[5] = t14;
    cgret[6] = (-t10 + t41 + t43 - t35 + t24 - t21 - t26 + t27) * t30;
    cgret[7] = (-t7 + t10 + t16 - t43 + t27 - t28 - t21 + t25) * t30;
    //cgret[8] = 1;
}

/// Compute homogaphy from 4 corresponding points.
void homography_from_4corresp(
    const double *a, const double *b, const double *c, const double *d,
    const double *x, const double *y, const double *z, const double *w,
    double R[3][3]) {
    double Hr[3][3], Hl[3][3];

    homography_from_4pt(a,b,c,d,&Hr[0][0]);
    homography_from_4pt(x,y,z,w,&Hl[0][0]);

    // the following code computes R = Hl * inverse Hr
    double t2 = Hr[1][1]-Hr[2][1]*Hr[1][2];
    double t4 = Hr[0][0]*Hr[1][1];
    double t5 = Hr[0][0]*Hr[1][2];
    double t7 = Hr[1][0]*Hr[0][1];
    double t8 = Hr[0][2]*Hr[1][0];
    double t10 = Hr[0][1]*Hr[2][0];
    double t12 = Hr[0][2]*Hr[2][0];
    double t15 = 1/(t4-t5*Hr[2][1]-t7+t8*Hr[2][1]+t10*Hr[1][2]-t12*Hr[1][1]);
    double t18 = -Hr[1][0]+Hr[1][2]*Hr[2][0];
    double t23 = -Hr[1][0]*Hr[2][1]+Hr[1][1]*Hr[2][0];
    double t28 = -Hr[0][1]+Hr[0][2]*Hr[2][1];
    double t31 = Hr[0][0]-t12;
    double t35 = Hr[0][0]*Hr[2][1]-t10;
    double t41 = -Hr[0][1]*Hr[1][2]+Hr[0][2]*Hr[1][1];
    double t44 = t5-t8;
    double t47 = t4-t7;
    double t48 = t2*t15;
    double t49 = t28*t15;
    double t50 = t41*t15;
    R[0][0] = Hl[0][0]*t48+Hl[0][1]*(t18*t15)-Hl[0][2]*(t23*t15);
    R[0][1] = Hl[0][0]*t49+Hl[0][1]*(t31*t15)-Hl[0][2]*(t35*t15);
    R[0][2] = -Hl[0][0]*t50-Hl[0][1]*(t44*t15)+Hl[0][2]*(t47*t15);
    R[1][0] = Hl[1][0]*t48+Hl[1][1]*(t18*t15)-Hl[1][2]*(t23*t15);
    R[1][1] = Hl[1][0]*t49+Hl[1][1]*(t31*t15)-Hl[1][2]*(t35*t15);
    R[1][2] = -Hl[1][0]*t50-Hl[1][1]*(t44*t15)+Hl[1][2]*(t47*t15);
    R[2][0] = Hl[2][0]*t48+Hl[2][1]*(t18*t15)-t23*t15;
    R[2][1] = Hl[2][0]*t49+Hl[2][1]*(t31*t15)-t35*t15;
    R[2][2] = -Hl[2][0]*t50-Hl[2][1]*(t44*t15)+t47*t15;
}
