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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#ifndef BERNOULLI_H
#define BERNOULLI_H

#include <math.h>
#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpz_vec.h"
#include "flint/fmpq.h"
#include "flint/arith.h"
#include "fmprb.h"
#include "arb.h"

#ifdef __cplusplus
extern "C" {
#endif

extern slong TLS_PREFIX bernoulli_cache_num;

extern TLS_PREFIX fmpq * bernoulli_cache;

void bernoulli_cache_compute(slong n);

/*
Crude bound for the bits in d(n) = denom(B_n).
By von Staudt-Clausen, d(n) = prod_{p-1 | n} p
    <= prod_{k | n} 2k
    <= n^{sigma_0(n)}.

We get a more accurate estimate taking the square root of this.
Further, at least for sufficiently large n,
sigma_0(n) < exp(1.066 log(n) / log(log(n))).
*/
static __inline__ slong bernoulli_denom_size(slong n)
{
    return 0.5 * 1.4427 * log(n) * pow(n, 1.066 / log(log(n)));
}

static __inline__ slong bernoulli_zeta_terms(ulong s, slong prec)
{
    slong N;
    N = pow(2.0, (prec + 1.0) / (s - 1.0));
    N += ((N % 2) == 0);
    return N;
}

static __inline__ slong bernoulli_power_prec(slong i, ulong s1, slong wp)
{
    slong p = wp - s1 * log(i) * 1.44269504088896341;
    return FLINT_MAX(p, 10);
}

/* we should technically add O(log(n)) guard bits, but this is unnecessary
   in practice since the denominator estimate is quite a bit larger
   than the true denominators
 */
static __inline__ slong bernoulli_global_prec(ulong nmax)
{
    return arith_bernoulli_number_size(nmax) + bernoulli_denom_size(nmax);
}


/* avoid potential numerical problems for very small n */
#define BERNOULLI_REV_MIN 32

typedef struct
{
    slong alloc;
    slong prec;
    slong max_power;
    fmpz * powers;
    fmpz_t pow_error;
    arb_t prefactor;
    arb_t two_pi_squared;
    ulong n;
}
bernoulli_rev_struct;

typedef bernoulli_rev_struct bernoulli_rev_t[1];

void bernoulli_rev_init(bernoulli_rev_t iter, ulong nmax);

void bernoulli_rev_next(fmpz_t numer, fmpz_t denom, bernoulli_rev_t iter);

void bernoulli_rev_clear(bernoulli_rev_t iter);


#define BERNOULLI_ENSURE_CACHED(n) \
  do { \
    slong __n = (n); \
    if (__n >= bernoulli_cache_num) \
        bernoulli_cache_compute(__n + 1); \
  } while (0); \

slong bernoulli_bound_2exp_si(ulong n);


void bernoulli_fmprb_ui_zeta(fmprb_t b, ulong n, slong prec);

void bernoulli_fmprb_ui(fmprb_t b, ulong n, slong prec);

void _bernoulli_fmpq_ui_zeta(fmpz_t num, fmpz_t den, ulong n);

void _bernoulli_fmpq_ui(fmpz_t num, fmpz_t den, ulong n);

void bernoulli_fmpq_ui(fmpq_t b, ulong n);

#ifdef __cplusplus
}
#endif

#endif

