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
#include "fmpcb.h"

extern TLS_PREFIX long * gamma_taylor_bound_mag_cache;
extern TLS_PREFIX fmpr_struct * gamma_taylor_bound_ratio_cache;
extern TLS_PREFIX long gamma_taylor_bound_cache_num;

extern TLS_PREFIX fmprb_ptr gamma_taylor_coeffs;
extern TLS_PREFIX long gamma_taylor_prec;
extern TLS_PREFIX long gamma_taylor_num;

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


void gamma_stirling_choose_param_fmprb(int * reflect, long * r, long * n, const fmprb_t z, int use_reflect, int digamma, long prec);
void gamma_stirling_choose_param_fmpcb(int * reflect, long * r, long * n, const fmpcb_t z, int use_reflect, int digamma, long prec);

void gamma_stirling_bound_phase(fmpr_t bound, const fmpcb_t z, long prec);
void gamma_stirling_bound_fmprb(fmpr_struct * err, const fmprb_t x, long k0, long knum, long n);
void gamma_stirling_bound_fmpcb(fmpr_struct * err, const fmpcb_t z, long k0, long knum, long n);

void gamma_stirling_coeff(fmprb_t b, ulong k, int digamma, long prec);
void gamma_stirling_eval_series_fmprb(fmprb_t s, const fmprb_t z, long nterms, int digamma, long prec);
void gamma_stirling_eval_series_fmpcb(fmpcb_t s, const fmpcb_t z, long nterms, int digamma, long prec);

void gamma_stirling_eval_fmprb_series(fmprb_ptr res, const fmprb_t z, long n, long num, long prec);
void gamma_stirling_eval_fmpcb_series(fmpcb_ptr res, const fmpcb_t z, long n, long num, long prec);

void gamma_rising_fmprb_ui_bsplit_simple(fmprb_t y, const fmprb_t x, ulong n, long prec);
void gamma_rising_fmprb_ui_bsplit_eight(fmprb_t y, const fmprb_t x, ulong n, long prec);
void gamma_rising_fmprb_ui_bsplit_rectangular(fmprb_t y, const fmprb_t x, ulong n, ulong step, long prec);
void gamma_rising_fmprb_ui_bsplit(fmprb_t y, const fmprb_t x, ulong n, long prec);

void gamma_rising_fmprb_ui_delta(fmprb_t y, const fmprb_t x, ulong n, ulong m, long prec);
void gamma_rising_fmpcb_ui_delta(fmpcb_t y, const fmpcb_t x, ulong n, ulong m, long prec);

void gamma_rising_fmpcb_ui_bsplit_simple(fmpcb_t y, const fmpcb_t x, ulong n, long prec);
void gamma_rising_fmpcb_ui_bsplit_eight(fmpcb_t y, const fmpcb_t x, ulong n, long prec);
void gamma_rising_fmpcb_ui_bsplit_rectangular(fmpcb_t y, const fmpcb_t x, ulong n, ulong step, long prec);
void gamma_rising_fmpcb_ui_bsplit(fmpcb_t y, const fmpcb_t x, ulong n, long prec);

void gamma_rising_fmprb_ui_multipoint(fmprb_t f, const fmprb_t c, ulong n, long prec);

void gamma_rising_fmprb_fmpq_ui_bsplit(fmprb_t y, const fmpq_t x, ulong n, long prec);

void gamma_harmonic_sum_fmprb_ui_bsplit_rectangular(fmprb_t y, const fmprb_t x, ulong n, ulong step, long prec);
void gamma_harmonic_sum_fmprb_ui_bsplit_simple(fmprb_t y, const fmprb_t x, ulong n, long prec);
void gamma_harmonic_sum_fmprb_ui_bsplit(fmprb_t y, const fmprb_t x, ulong n, long prec);

void gamma_harmonic_sum_fmpcb_ui_bsplit_rectangular(fmpcb_t y, const fmpcb_t x, ulong n, ulong step, long prec);
void gamma_harmonic_sum_fmpcb_ui_bsplit_simple(fmpcb_t y, const fmpcb_t x, ulong n, long prec);
void gamma_harmonic_sum_fmpcb_ui_bsplit(fmpcb_t y, const fmpcb_t x, ulong n, long prec);

void gamma_lgamma_series_at_one(fmprb_ptr u, long len, long prec);
void gamma_series_fmpq_hypgeom(fmprb_ptr res, const fmpq_t a, long len, long prec);
void gamma_small_frac(fmprb_t y, unsigned int p, unsigned int q, long prec);

void gamma_fmpq_outward(fmprb_t y, const fmpq_t x, long prec);
void gamma_fmpq_stirling(fmprb_t y, const fmpq_t a, long prec);

#endif

