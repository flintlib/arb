/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef ARB_FPWRAP_H
#define ARB_FPWRAP_H

#ifdef ARB_FPWRAP_INLINES_C
#define ARB_FPWRAP_INLINE
#else
#define ARB_FPWRAP_INLINE static __inline__
#endif

#include "arb.h"
#include "acb.h"

#ifdef __cplusplus
extern "C" {
#endif

#define FPWRAP_SUCCESS 0
#define FPWRAP_UNABLE 1

#define FPWRAP_ACCURATE_PARTS 1
#define FPWRAP_CORRECT_ROUNDING 2
#define FPWRAP_WORK_LIMIT 65536

typedef struct
{
    double real;
    double imag;
}
complex_double;

int arb_fpwrap_double_exp(double * res, double x, int flags);
int arb_fpwrap_cdouble_exp(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_expm1(double * res, double x, int flags);
int arb_fpwrap_cdouble_expm1(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_log(double * res, double x, int flags);
int arb_fpwrap_cdouble_log(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_log1p(double * res, double x, int flags);
int arb_fpwrap_cdouble_log1p(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_sqrt(double * res, double x, int flags);
int arb_fpwrap_cdouble_sqrt(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_rsqrt(double * res, double x, int flags);
int arb_fpwrap_cdouble_rsqrt(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_cbrt(double * res, double x, int flags);
int arb_fpwrap_cdouble_cbrt(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_sin(double * res, double x, int flags);
int arb_fpwrap_cdouble_sin(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_cos(double * res, double x, int flags);
int arb_fpwrap_cdouble_cos(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_tan(double * res, double x, int flags);
int arb_fpwrap_cdouble_tan(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_cot(double * res, double x, int flags);
int arb_fpwrap_cdouble_cot(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_sec(double * res, double x, int flags);
int arb_fpwrap_cdouble_sec(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_csc(double * res, double x, int flags);
int arb_fpwrap_cdouble_csc(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_sinc(double * res, double x, int flags);
int arb_fpwrap_cdouble_sinc(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_sin_pi(double * res, double x, int flags);
int arb_fpwrap_cdouble_sin_pi(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_cos_pi(double * res, double x, int flags);
int arb_fpwrap_cdouble_cos_pi(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_tan_pi(double * res, double x, int flags);
int arb_fpwrap_cdouble_tan_pi(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_cot_pi(double * res, double x, int flags);
int arb_fpwrap_cdouble_cot_pi(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_sinc_pi(double * res, double x, int flags);
int arb_fpwrap_cdouble_sinc_pi(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_asin(double * res, double x, int flags);
int arb_fpwrap_cdouble_asin(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_acos(double * res, double x, int flags);
int arb_fpwrap_cdouble_acos(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_atan(double * res, double x, int flags);
int arb_fpwrap_cdouble_atan(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_atan2(double * res, double x1, double x2, int flags);

int arb_fpwrap_double_asinh(double * res, double x, int flags);
int arb_fpwrap_cdouble_asinh(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_acosh(double * res, double x, int flags);
int arb_fpwrap_cdouble_acosh(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_atanh(double * res, double x, int flags);
int arb_fpwrap_cdouble_atanh(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_rising(double * res, double x, double n, int flags);
int arb_fpwrap_cdouble_rising(complex_double * res, complex_double x, complex_double n, int flags);

int arb_fpwrap_double_gamma(double * res, double x, int flags);
int arb_fpwrap_cdouble_gamma(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_rgamma(double * res, double x, int flags);
int arb_fpwrap_cdouble_rgamma(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_lgamma(double * res, double x, int flags);
int arb_fpwrap_cdouble_lgamma(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_digamma(double * res, double x, int flags);
int arb_fpwrap_cdouble_digamma(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_zeta(double * res, double x, int flags);
int arb_fpwrap_cdouble_zeta(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_hurwitz_zeta(double * res, double s, double z, int flags);
int arb_fpwrap_cdouble_hurwitz_zeta(complex_double * res, complex_double s, complex_double z, int flags);

int arb_fpwrap_double_barnes_g(double * res, double x, int flags);
int arb_fpwrap_cdouble_barnes_g(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_log_barnes_g(double * res, double x, int flags);
int arb_fpwrap_cdouble_log_barnes_g(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_polygamma(double * res, double s, double z, int flags);
int arb_fpwrap_cdouble_polygamma(complex_double * res, complex_double s, complex_double z, int flags);

int arb_fpwrap_double_polylog(double * res, double s, double z, int flags);
int arb_fpwrap_cdouble_polylog(complex_double * res, complex_double s, complex_double z, int flags);

int arb_fpwrap_double_agm1(double * res, double x, int flags);
int arb_fpwrap_cdouble_agm1(complex_double * res, complex_double x, int flags);

int arb_fpwrap_double_agm(double * res, double x, double y, int flags);
int arb_fpwrap_cdouble_agm(complex_double * res, complex_double x, complex_double y, int flags);



#ifdef __cplusplus
}
#endif

#endif

