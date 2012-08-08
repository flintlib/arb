/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#ifndef UFLOAT_H
#define UFLOAT_H

#include <stdio.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"
#include "math.h"

static __inline__ mp_limb_t
n_rshift_ceil(mp_limb_t a, int k)
{
    return (a >> k) + (((a >> k) << k) != a);
}

/*
A ufloat holds a nonzero bound for the magnitude of a number.
The mantissa has exactly UFLOAT_PREC bits, which must be at most
FLINT_BITS / 2.
*/

#define UFLOAT_PREC 15
#define UFLOAT_MAN_MIN (1UL << (UFLOAT_PREC - 1))
#define UFLOAT_MAN_MAX ((1UL << UFLOAT_PREC) - 1)

typedef struct
{
    mp_limb_t man;
    long exp;
}
ufloat_struct;

typedef ufloat_struct ufloat_t[1];

static __inline__ void
ufloat_print(const ufloat_t x)
{
    double t = ldexp(x->man, x->exp - UFLOAT_PREC);

    printf("ufloat{man=%lu, exp=%ld, %.20g}\n", x->man, x->exp, t);
}

static __inline__ void
ufloat_normalise(ufloat_t z)
{
    int b = FLINT_BIT_COUNT(z->man);

    if (b < UFLOAT_PREC)
    {
        z->man <<= (UFLOAT_PREC - b);
        z->exp -= (UFLOAT_PREC - b);
    }
    else
    {
        int adjust;

        z->man = n_rshift_ceil(z->man, (b - UFLOAT_PREC));
        z->exp += (b - UFLOAT_PREC);

        /* possible overflow to power of two */
        adjust = z->man >> UFLOAT_PREC;
        z->man = n_rshift_ceil(z->man, adjust);
        z->exp += adjust;
    }
}


static __inline__ int
ufloat_cmp_one(const ufloat_t u)
{
    if (u->exp == 1L)
        return u->man == (1UL << (UFLOAT_PREC - 1)) ? 0 : 1;
    else
        return u->exp < 1L ? -1 : 1;
}

static __inline__ void
ufloat_set(ufloat_t u, const ufloat_t v)
{
    u->man = v->man;
    u->exp = v->exp;
}


static __inline__ void
ufloat_zero(ufloat_t u)
{
    u->man = 0UL;
    u->exp = 0L;
}

static __inline__ int
ufloat_is_zero(ufloat_t u)
{
    return u->man != 0UL;
}


/* TODO: add constant_p handling; overflow */
static __inline__ void
ufloat_set_ui_2exp(ufloat_t u, mp_limb_t a, long exp)
{
    u->man = a;
    u->exp = exp + UFLOAT_PREC;
    ufloat_normalise(u);
}

/* TODO: overflow */
static __inline__ void
ufloat_set_ll_2exp(ufloat_t u, mp_limb_t hi, mp_limb_t lo, long exp)
{
    if (hi == 0UL)
    {
        ufloat_set_ui_2exp(u, lo, exp);
    }
    else
    {
        unsigned int c;
        count_leading_zeros(c, hi);

        if (c != 0)
            hi = (hi << c) | (lo >> (FLINT_BITS - c));

        if (hi + 1UL != 0UL)
            ufloat_set_ui_2exp(u, hi + 1UL, exp + FLINT_BITS - c);
        else
            ufloat_set_ui_2exp(u, 1UL, exp + 2 * FLINT_BITS - c);
    }
}

static __inline__ void
ufloat_one(ufloat_t u)
{
    u->man = (1UL << (UFLOAT_PREC - 1));
    u->exp = 1L;
}

static __inline__ void
ufloat_add(ufloat_t z, const ufloat_t x, const ufloat_t y)
{
    long shift;
    int adjust;

    shift = x->exp - y->exp;

    if (shift >= 0)
    {
        z->exp = x->exp;
        if (shift >= UFLOAT_PREC)
            z->man = x->man + 1;
        else
            z->man = x->man + n_rshift_ceil(y->man, shift);
    }
    else
    {
        shift = -shift;
        z->exp = y->exp;
        if (shift >= UFLOAT_PREC)
            z->man = y->man + 1;
        else
            z->man = y->man + n_rshift_ceil(x->man, shift);
    }

    /* adjust for carry */
    adjust = z->man >> UFLOAT_PREC;
    z->man = n_rshift_ceil(z->man, adjust);
    z->exp += adjust;
}

/* bound |x-y| */
static __inline__ void
ufloat_sub(ufloat_t z, const ufloat_t x, const ufloat_t y)
{
    long shift;

    mp_limb_t a, b;

    shift = x->exp - y->exp;

    if (shift == 0)
    {
        if (x->man >= y->man)
            z->man = x->man - y->man;
        else
            z->man = y->man - x->man;
        z->exp = x->exp;
        ufloat_normalise(z);
        return;
    }

    if (shift > 0)
    {
        z->exp = x->exp;
        a = x->man;
        b = y->man;
    }
    else
    {
        shift = -shift;
        z->exp = y->exp;
        a = y->man;
        b = x->man;
    }

    if (shift >= UFLOAT_PREC)
    {
        z->man = a;
    }
    else
    {
        z->man = (a << shift) - b;
        z->exp -= shift;
    }

    ufloat_normalise(z);
}

static __inline__ void
ufloat_set_2exp(ufloat_t y, const ufloat_t x, long exp)
{

}

static __inline__ void
ufloat_add_2exp(ufloat_t y, const ufloat_t x, long exp)
{
    long s, shift = x->exp - exp;

    /* x is larger than 2^exp */
    if (shift > 0)
    {
        s = FLINT_MAX(0L, UFLOAT_PREC - shift);
        y->man = x->man + (1UL << s);
        y->exp = x->exp;
    }
    /* 2^exp is larger than x */
    else
    {
        s = FLINT_MIN(UFLOAT_PREC - 1L, 1L - shift);
        y->man = (1UL << (UFLOAT_PREC - 1)) + (x->man >> s) + 1UL;
        y->exp = exp + 1L;
    }

    /* unlikely */
    if (y->man > UFLOAT_MAN_MAX)
    {
        y->man = n_rshift_ceil(y->man, 1);
        y->exp++;
    }
}


static __inline__ void
ufloat_mul(ufloat_t z, const ufloat_t x, const ufloat_t y)
{
    mp_limb_t m;
    int adjust, shift;

    z->exp = x->exp + y->exp;
    m = x->man * y->man;

    /* product has one bit too little */
    adjust = m >> (2 * UFLOAT_PREC - 1);
    shift = UFLOAT_PREC - 1 + adjust;
    z->exp += shift;

    m = n_rshift_ceil(m, shift);

    /* power of two */
    adjust = m >> UFLOAT_PREC;
    m >>= adjust;
    z->exp += adjust;

    z->man = m;
}

static __inline__ void
ufloat_addmul(ufloat_t z, const ufloat_t x, const ufloat_t y)
{
    ufloat_t t;
    ufloat_mul(t, x, y);
    ufloat_add(z, z, t);
}

static __inline__ void
ufloat_div(ufloat_t z, const ufloat_t x, const ufloat_t y)
{
    int adjust;

    z->man = (x->man << UFLOAT_PREC) / y->man + 1;
    z->exp = x->exp - y->exp;

    /* quotient can be one bit too large */
    adjust = z->man >> UFLOAT_PREC;
    z->man = n_rshift_ceil(z->man, adjust);
    z->exp += adjust;

    /* power of two */
    adjust = z->man >> UFLOAT_PREC;
    z->man >>= adjust;
    z->exp += adjust;
}

static __inline__ void
ufloat_randtest(ufloat_t u, flint_rand_t state, long erange)
{
    u->man = n_randbits(state, UFLOAT_PREC);
    u->man |= (1UL << (UFLOAT_PREC - 1));
    u->exp = n_randint(state, erange) - (erange / 2) + 1;
}

static __inline__ void
ufloat_set_fmpz(ufloat_t u, const fmpz_t z)
{
    u->man = fmpz_abs_ubound_ui_2exp(&u->exp, z, UFLOAT_PREC);
    u->exp += UFLOAT_PREC;
}

static __inline__ void
ufloat_set_fmpz_lower(ufloat_t u, const fmpz_t z)
{
    u->man = fmpz_abs_lbound_ui_2exp(&u->exp, z, UFLOAT_PREC);
    u->exp += UFLOAT_PREC;
}

static __inline__ void
ufloat_get_fmpz(fmpz_t z, const ufloat_t u)
{
    long exp = u->exp - UFLOAT_PREC;

    _fmpz_demote(z);
    *z = u->man;

    if (exp >= 0)
    {
        if (exp < FLINT_BITS - UFLOAT_PREC - 2)
            *z = (u->man << exp);
        else
            fmpz_mul_2exp(z, z, exp);
    }
    else
    {
        if (exp > -FLINT_BITS)
            *z = n_rshift_ceil(u->man, -exp);
        else
            *z = 1;
    }
}

static __inline__ int
ufloat_equal(const ufloat_t u, const ufloat_t v)
{
    return u->man == v->man && u->exp == v->exp;
}

void ufloat_get_mpfr(mpfr_t x, const ufloat_t u);

void ufloat_set_mpfr(ufloat_t u, const mpfr_t x);

void ufloat_log(ufloat_t z, const ufloat_t x);

void ufloat_log1p(ufloat_t z, const ufloat_t x);

#endif
