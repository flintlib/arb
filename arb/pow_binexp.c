/*
    Copyright (C) 2012, 2013 Fredrik Johansson
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

void
_arb_pow_mpz_binexp(arb_t y, const arb_t b, mpz_srcptr e, slong prec)
{
    slong i, wp, bits;

    if (e->_mp_size < 0)
    {
        __mpz_struct f = { 1, 0, NULL };
        f._mp_size = -(e->_mp_size);
        f._mp_d = e->_mp_d;

        if (arb_is_exact(b))
        {
            _arb_pow_mpz_binexp(y, b, &f, prec + 2);
            arb_inv(y, y, prec);
        }
        else
        {
            arb_inv(y, b, prec + mpz_sizeinbase(e, 2) + 2);
            _arb_pow_mpz_binexp(y, y, &f, prec);
        }

        return;
    }

    if (y == b)
    {
        arb_t t;
        arb_init(t);
        arb_set(t, b);
        _arb_pow_mpz_binexp(y, t, e, prec);
        arb_clear(t);
        return;
    }

    arb_set(y, b);

    bits = mpz_sizeinbase(e, 2);
    wp = ARF_PREC_ADD(prec, bits);

    for (i = bits - 2; i >= 0; i--)
    {
        arb_sqr(y, y, wp);
        if (mpz_tstbit(e, i))
            arb_mul(y, y, b, wp);
    }
}

void
arb_pow_ui_binexp(arb_t y, const arb_t b, ulong e, slong prec)
{
    slong i, wp, bits;

    if (e <= 2)
    {
        if (e == 0)
            arb_one(y);
        else if (e == 1)
            arb_set_round(y, b, prec);
        else
            arb_sqr(y, b, prec);
        return;
    }

    if (y == b)
    {
        arb_t t;
        arb_init(t);
        arb_set(t, b);
        arb_pow_ui_binexp(y, t, e, prec);
        arb_clear(t);
        return;
    }

    arb_set(y, b);

    bits = FLINT_BIT_COUNT(e);
    wp = ARF_PREC_ADD(prec, bits);

    for (i = bits - 2; i >= 0; i--)
    {
        arb_sqr(y, y, wp);
        if (((WORD(1) << i) & e) != 0)
            arb_mul(y, y, b, wp);
    }
}
