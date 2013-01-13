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
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq.h"
#include "arith.h"
#include "fmprb.h"

extern long __thread bernoulli_cache_num;

extern __thread fmpq * bernoulli_cache;

void bernoulli_cache_compute(long n);

/*
Crude bound for the bits in d(n) = denom(B_n).
By von Staudt-Clausen, d(n) = prod_{p-1 | n} p
    <= prod_{k | n} 2k
    <= n^{sigma_0(n)}.

We get a more accurate estimate taking the square root of this.
Further, at least for sufficiently large n,
sigma_0(n) < exp(1.066 log(n) / log(log(n))).
*/
static __inline__ long denom_size(long n)
{
    return 0.5 * 1.4427 * log(n) * pow(n, 1.066 / log(log(n)));
}

static __inline__ long zeta_terms(ulong s, long prec)
{
    long N;
    N = pow(2.0, (prec + 1.0) / (s - 1.0));
    N += ((N % 2) == 0);
    return N;
}

static __inline__ long power_prec(long i, ulong s1, long wp)
{
    long p = wp - s1 * log(i) * 1.44269504088896341;
    return FLINT_MAX(p, 10);
}

/* we should technically add O(log(n)) guard bits, but this is unnecessary
   in practice since the denominator estimate is quite a bit larger
   than the true denominators
 */
static __inline__ long global_prec(ulong nmax)
{
    return arith_bernoulli_number_size(nmax) + denom_size(nmax);
}


/* avoid potential numerical problems for very small n */
#define bernoulli_rev_MIN 32

typedef struct
{
    long alloc;
    long prec;
    long max_power;
    fmpz * powers;
    fmpz_t pow_error;
    fmprb_t prefactor;
    fmprb_t two_pi_squared;
    ulong n;
}
bernoulli_rev_struct;

typedef bernoulli_rev_struct bernoulli_rev_t[1];

void bernoulli_rev_init(bernoulli_rev_t iter, ulong nmax);

void bernoulli_rev_next(fmpz_t numer, fmpz_t denom, bernoulli_rev_t iter);

void bernoulli_rev_clear(bernoulli_rev_t iter);

#endif

