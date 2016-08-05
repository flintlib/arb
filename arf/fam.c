/*
    Copyright (C) 2014 Fredrik Johansson
    Copyright (C) 2016 Ricky E. Farr

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arf.h"

int
arf_fam(arf_t w, const arf_t z, const arf_t x, const arf_t y,
        slong prec, arf_rnd_t rnd)
{
    if (arf_is_special(x) || arf_is_special(y) || arf_is_special(z))
    {
        if (arf_is_zero(z))
        {
            return arf_mul(w, x, y, prec, rnd);
        }
        else if (arf_is_finite(x) && arf_is_finite(y))
        {
            return arf_set_round(w, z, prec, rnd);
        }
        else
        {
            if (w != z)
            {
                arf_mul(w, x, y, ARF_PREC_EXACT, ARF_RND_DOWN);
                return arf_add(w, z, w, prec, rnd);
            }

            /* todo: speed up */
            int inexact;
            arf_t t;
            arf_init(t);
            arf_mul(t, x, y, ARF_PREC_EXACT, ARF_RND_DOWN);
            inexact = arf_add(w, z, t, prec, rnd);
            arf_clear(t);
            return inexact;
        }
    }

    mp_size_t xn, yn, zn, tn, alloc;
    mp_srcptr xptr, yptr, zptr;
    mp_ptr tptr;
    fmpz_t texp;
    slong shift;
    int tsgnbit, inexact;
    ARF_MUL_TMP_DECL

    tsgnbit = ARF_SGNBIT(x) ^ ARF_SGNBIT(y);
    ARF_GET_MPN_READONLY(xptr, xn, x);
    ARF_GET_MPN_READONLY(yptr, yn, y);
    ARF_GET_MPN_READONLY(zptr, zn, z);

    fmpz_init(texp);

    _fmpz_add2_fast(texp, ARF_EXPREF(x), ARF_EXPREF(y), 0);

    shift = _fmpz_sub_small(ARF_EXPREF(z), texp);

    alloc = tn = xn + yn;
    ARF_MUL_TMP_ALLOC(tptr, alloc);

    ARF_MPN_MUL(tptr, xptr, xn, yptr, yn);

    tn   -= (tptr[0] == 0);
    tptr += (tptr[0] == 0);

    if (shift >= 0)
        inexact = _arf_add_mpn(w, zptr, zn, ARF_SGNBIT(z), ARF_EXPREF(z),
            tptr, tn, tsgnbit, shift, prec, rnd);
    else
        inexact = _arf_add_mpn(w, tptr, tn, tsgnbit, texp,
            zptr, zn, ARF_SGNBIT(z), -shift, prec, rnd);

    ARF_MUL_TMP_FREE(tptr, alloc)
    fmpz_clear(texp);

    return inexact;
}

int
arf_fam_mpz(arf_ptr w, arf_srcptr z, arf_srcptr x, const mpz_t y,
            slong prec, arf_rnd_t rnd)
{
    mp_size_t yn = FLINT_ABS(y->_mp_size);

    if (arf_is_special(x) || yn == 0 || arf_is_special(z))
    {
        if (arf_is_zero(z))
        {
            return arf_mul_mpz(w, x, y, prec, rnd);
        }
        else if (arf_is_finite(x))
        {
            return arf_set_round(w, z, prec, rnd);
        }
        else
        {
            /* todo: speed up */
            arf_mul_mpz(w, x, y, ARF_PREC_EXACT, ARF_RND_DOWN);
            return arf_add(w, z, w, prec, rnd);
        }
    }

    mp_size_t xn, zn, tn, alloc;
    mp_srcptr xptr, yptr, zptr;
    mp_ptr tptr;
    fmpz_t texp, yexp;
    slong shift;
    int tsgnbit, ysgnbit, inexact;
    ARF_MUL_TMP_DECL

    ARF_GET_MPN_READONLY(xptr, xn, x);

    yptr = y->_mp_d;
    ysgnbit = (y->_mp_size < 0);
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
        inexact = _arf_add_mpn(w, zptr, zn, ARF_SGNBIT(z), ARF_EXPREF(z),
            tptr, tn, tsgnbit, shift, prec, rnd);
    else
        inexact = _arf_add_mpn(w, tptr, tn, tsgnbit, texp,
            zptr, zn, ARF_SGNBIT(z), -shift, prec, rnd);

    ARF_MUL_TMP_FREE(tptr, alloc)
    fmpz_clear(texp);

    return inexact;
}
