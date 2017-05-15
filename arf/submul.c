/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arf.h"

int
arf_submul(arf_ptr z, arf_srcptr x, arf_srcptr y, slong prec, arf_rnd_t rnd)
{
    mp_size_t xn, yn, zn, tn, alloc;
    mp_srcptr xptr, yptr, zptr;
    mp_ptr tptr;
    fmpz_t texp;
    slong shift;
    int tsgnbit, inexact;
    ARF_MUL_TMP_DECL

    if (arf_is_special(x) || arf_is_special(y) || arf_is_special(z))
    {
        if (arf_is_zero(z))
        {
            return arf_neg_mul(z, x, y, prec, rnd);
        }
        else if (arf_is_finite(x) && arf_is_finite(y))
        {
            return arf_set_round(z, z, prec, rnd);
        }
        else
        {
            /* todo: speed up */
            arf_t t;
            arf_init(t);
            arf_mul(t, x, y, ARF_PREC_EXACT, ARF_RND_DOWN);
            inexact = arf_sub(z, z, t, prec, rnd);
            arf_clear(t);
            return inexact;
        }
    }

    tsgnbit = ARF_SGNBIT(x) ^ ARF_SGNBIT(y) ^ 1;
    ARF_GET_MPN_READONLY(xptr, xn, x);
    ARF_GET_MPN_READONLY(yptr, yn, y);
    ARF_GET_MPN_READONLY(zptr, zn, z);

    fmpz_init(texp);

    _fmpz_add2_fast(texp, ARF_EXPREF(x), ARF_EXPREF(y), 0);
    shift = _fmpz_sub_small(ARF_EXPREF(z), texp);

    alloc = tn = xn + yn;
    ARF_MUL_TMP_ALLOC(tptr, alloc)

    ARF_MPN_MUL(tptr, xptr, xn, yptr, yn);

    tn -= (tptr[0] == 0);
    tptr += (tptr[0] == 0);

    if (shift >= 0)
        inexact = _arf_add_mpn(z, zptr, zn, ARF_SGNBIT(z), ARF_EXPREF(z),
            tptr, tn, tsgnbit, shift, prec, rnd);
    else
        inexact = _arf_add_mpn(z, tptr, tn, tsgnbit, texp,
            zptr, zn, ARF_SGNBIT(z), -shift, prec, rnd);

    ARF_MUL_TMP_FREE(tptr, alloc)
    fmpz_clear(texp);

    return inexact;
}

int
arf_submul_mpz(arf_ptr z, arf_srcptr x, const mpz_t y, slong prec, arf_rnd_t rnd)
{
    mp_size_t xn, yn, zn, tn, alloc;
    mp_srcptr xptr, yptr, zptr;
    mp_ptr tptr;
    fmpz_t texp, yexp;
    slong shift;
    int tsgnbit, ysgnbit, inexact;
    ARF_MUL_TMP_DECL

    yn = FLINT_ABS(y->_mp_size);

    if (arf_is_special(x) || yn == 0 || arf_is_special(z))
    {
        if (arf_is_zero(z))
        {
            /* TODO: make more efficient */
            arf_mul_mpz(z, x, y, ARF_PREC_EXACT, rnd);
            return arf_neg_round(z, z, prec, rnd);
        }
        else if (arf_is_finite(x))
        {
            return arf_set_round(z, z, prec, rnd);
        }
        else
        {
            /* todo: speed up */
            arf_t t;
            arf_init(t);
            arf_mul_mpz(t, x, y, ARF_PREC_EXACT, ARF_RND_DOWN);
            inexact = arf_sub(z, z, t, prec, rnd);
            arf_clear(t);
            return inexact;
        }
    }

    ARF_GET_MPN_READONLY(xptr, xn, x);

    yptr = y->_mp_d;
    ysgnbit = (y->_mp_size > 0);
    *yexp = yn * FLINT_BITS;

    ARF_GET_MPN_READONLY(zptr, zn, z);

    fmpz_init(texp);

    tsgnbit = ARF_SGNBIT(x) ^ ysgnbit;

    alloc = tn = xn + yn;
    ARF_MUL_TMP_ALLOC(tptr, alloc)

    ARF_MPN_MUL(tptr, xptr, xn, yptr, yn);

    shift = (tptr[tn - 1] == 0) * FLINT_BITS;
    tn -= (tptr[tn - 1] == 0);

    _fmpz_add2_fast(texp, ARF_EXPREF(x), yexp, -shift);
    shift = _fmpz_sub_small(ARF_EXPREF(z), texp);

    if (shift >= 0)
        inexact = _arf_add_mpn(z, zptr, zn, ARF_SGNBIT(z), ARF_EXPREF(z),
            tptr, tn, tsgnbit, shift, prec, rnd);
    else
        inexact = _arf_add_mpn(z, tptr, tn, tsgnbit, texp,
            zptr, zn, ARF_SGNBIT(z), -shift, prec, rnd);

    ARF_MUL_TMP_FREE(tptr, alloc)
    fmpz_clear(texp);

    return inexact;
}

