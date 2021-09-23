/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <math.h>

#include "arb.h"
#include "acb.h"
#include "arb_fpwrap.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"


int
arb_accurate_enough_d(const arb_t x, int flags)
{
    if (flags & FPWRAP_CORRECT_ROUNDING)
    {
        return arb_can_round_arf(x, 53, ARF_RND_NEAR);
    }

    if (arb_rel_accuracy_bits(x) >= 53 + 1)
        return 1;

    /* Rounding will give +/- 0 (we don't worry which) */
    if (mag_cmp_2exp_si(arb_radref(x), -1077) < 0 &&
            arf_cmpabs_2exp_si(arb_midref(x), -1077) < 0)
    {
        return 1;
    }

    /* Rounding will give +/- inf */
    if (arb_rel_accuracy_bits(x) > 2 &&
        arf_cmpabs_2exp_si(arb_midref(x), 1024) > 0)
    {
        return 1;
    }

    return 0;
}

int
acb_accurate_enough_d(const acb_t x, int flags)
{
    if (flags & FPWRAP_CORRECT_ROUNDING)
    {
        return arb_can_round_arf(acb_realref(x), 53, ARF_RND_NEAR) &&
               arb_can_round_arf(acb_imagref(x), 53, ARF_RND_NEAR);
    }

    if (flags & FPWRAP_ACCURATE_PARTS)
    {
        return arb_accurate_enough_d(acb_realref(x), flags) &&
               arb_accurate_enough_d(acb_imagref(x), flags);
    }

    if (acb_rel_accuracy_bits(x) >= 53 + 1)
        return 1;

    /* Rounding will give +/- 0 (we don't worry which) */
    if (mag_cmp_2exp_si(arb_radref(acb_realref(x)), -1077) < 0 &&
            arf_cmpabs_2exp_si(arb_midref(acb_realref(x)), -1077) < 0 &&
            mag_cmp_2exp_si(arb_radref(acb_imagref(x)), -1077) < 0 &&
            arf_cmpabs_2exp_si(arb_midref(acb_imagref(x)), -1077) < 0)
    {
        return 1;
    }

    /* Rounding will give +/- inf */
    if (acb_rel_accuracy_bits(x) > 2 &&
        (arf_cmpabs_2exp_si(arb_midref(acb_realref(x)), 1024) > 0 ||
         arf_cmpabs_2exp_si(arb_midref(acb_imagref(x)), 1024) > 0))
    {
        return 1;
    }

    return 0;
}

#define WP_INITIAL 64

static slong
double_wp_max(int flags)
{
    int iters;

    iters = flags / FPWRAP_WORK_LIMIT;

    if (iters <= 0)
        return 64 << 7;

    if (iters >= 25)
        return 64 << 24;

    return 64 << iters;
}


typedef void (*arb_func_1)(arb_t, const arb_t, slong prec);
typedef void (*acb_func_1)(acb_t, const acb_t, slong prec);

int arb_fpwrap_double(double * res, arb_func_1 func, double x, int flags)
{
    arb_t arb_res, arb_x;
    slong wp;
    int status;

    arb_init(arb_res);
    arb_init(arb_x);

    arb_set_d(arb_x, x);

    if (!arb_is_finite(arb_x))
    {
        arb_clear(arb_x);
        *res = D_NAN;
        return FPWRAP_UNABLE;
    }

    for (wp = WP_INITIAL; ; wp *= 2)
    {
        func(arb_res, arb_x, wp);

        if (arb_accurate_enough_d(arb_res, flags))
        {
            *res = arf_get_d(arb_midref(arb_res), ARF_RND_NEAR);
            status = FPWRAP_SUCCESS;
            break;
        }

        if (wp >= double_wp_max(flags))
        {
            *res = D_NAN;
            status = FPWRAP_UNABLE;
            break;
        }
    }

    arb_clear(arb_x);
    arb_clear(arb_res);

    return status;
}

int arb_fpwrap_cdouble(complex_double * res, acb_func_1 func, complex_double x, int flags)
{
    acb_t acb_res, acb_x;
    slong wp;
    int status;

    acb_init(acb_res);
    acb_init(acb_x);

    acb_set_d_d(acb_x, x.real, x.imag);

    /* no need to clear */
    if (!acb_is_finite(acb_x))
    {
        acb_clear(acb_x);
        res->real = D_NAN;
        res->imag = D_NAN;
        return FPWRAP_UNABLE;
    }

    for (wp = WP_INITIAL; ; wp *= 2)
    {
        func(acb_res, acb_x, wp);

        if (acb_accurate_enough_d(acb_res, flags))
        {
            res->real = arf_get_d(arb_midref(acb_realref(acb_res)), ARF_RND_NEAR);
            res->imag = arf_get_d(arb_midref(acb_imagref(acb_res)), ARF_RND_NEAR);
            status = FPWRAP_SUCCESS;
            break;
        }

        if (wp >= double_wp_max(flags))
        {
            res->real = D_NAN;
            res->imag = D_NAN;
            status = FPWRAP_UNABLE;
            break;
        }
    }

    acb_clear(acb_x);
    acb_clear(acb_res);

    return status;
}

#define DEF_DOUBLE_FUN(name, arb_fun) \
    int arb_fpwrap_double_ ## name(double * res, double x, int flags) \
    { \
        return arb_fpwrap_double(res, arb_fun, x, flags); \
    } \

#define DEF_CDOUBLE_FUN(name, acb_fun) \
    int arb_fpwrap_cdouble_ ## name(complex_double * res, complex_double x, int flags) \
    { \
        return arb_fpwrap_cdouble(res, acb_fun, x, flags); \
    } \


DEF_DOUBLE_FUN(gamma, arb_gamma)
DEF_CDOUBLE_FUN(gamma, acb_gamma)

DEF_DOUBLE_FUN(rgamma, arb_rgamma)
DEF_CDOUBLE_FUN(rgamma, acb_rgamma)

DEF_DOUBLE_FUN(lgamma, arb_lgamma)
DEF_CDOUBLE_FUN(lgamma, acb_lgamma)

DEF_DOUBLE_FUN(digamma, arb_digamma)
DEF_CDOUBLE_FUN(digamma, acb_digamma)

DEF_DOUBLE_FUN(zeta, arb_zeta)
DEF_CDOUBLE_FUN(zeta, acb_zeta)


