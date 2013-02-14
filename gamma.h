/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#ifndef GAMMA_H
#define GAMMA_H

#include <math.h>
#include "flint.h"
#include "fmprb.h"

extern __thread long * gamma_taylor_bound_mag_cache;
extern __thread fmpr_struct * gamma_taylor_bound_ratio_cache;
extern __thread long gamma_taylor_bound_cache_num;

extern __thread fmprb_struct * gamma_taylor_coeffs;
extern __thread long gamma_taylor_prec;
extern __thread long gamma_taylor_num;

void gamma_taylor_bound_ratio(fmpr_t r, long n);

long gamma_taylor_bound_mag(long n);

void gamma_taylor_bound_extend_cache(long n);

static __inline__ long
gamma_taylor_bound_cached(long n)
{
    if (n >= gamma_taylor_bound_cache_num)
        gamma_taylor_bound_extend_cache(n);

    return gamma_taylor_bound_mag_cache[n];
}

void gamma_taylor_bound_remainder(fmpr_t err, const fmpr_t z, long n);

static __inline__ long
gamma_taylor_coeffs_for_prec(long prec)
{
    long n = 5;
    while (gamma_taylor_bound_cached(n) - n > -prec)
        n++;
    return n;
}

void gamma_taylor_precompute(long num, long prec);

void gamma_taylor_eval_series_fmprb(fmprb_t y, const fmprb_t x, long prec);

void gamma_taylor_fmprb(fmprb_t y, const fmprb_t x, long prec);

#endif

