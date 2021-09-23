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


#ifdef __cplusplus
}
#endif

#endif

