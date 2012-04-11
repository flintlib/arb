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
    double t = ldexp(x->man, x->exp);

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
    z->exp = x->exp - UFLOAT_PREC - y->exp;

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
}

static __inline__ void
ufloat_set_fmpz_lower(ufloat_t u, const fmpz_t z)
{
    u->man = fmpz_abs_lbound_ui_2exp(&u->exp, z, UFLOAT_PREC);
}

static __inline__ void
ufloat_get_fmpz(fmpz_t z, const ufloat_t u)
{
    _fmpz_demote(z);
    *z = u->man;

    if (u->exp >= 0)
    {
        if (u->exp < FLINT_BITS - UFLOAT_PREC - 2)
            *z = (u->man << u->exp);
        else
            fmpz_mul_2exp(z, z, u->exp);
    }
    else
    {
        if (u->exp > -FLINT_BITS)
            *z = n_rshift_ceil(u->man, -u->exp);
        else
            *z = 1;
    }
}

void ufloat_log(ufloat_t z, const ufloat_t x);

void ufloat_log1p(ufloat_t z, const ufloat_t x);

#endif
