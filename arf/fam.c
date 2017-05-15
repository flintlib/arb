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
    if (w == z)
    {
        return arf_addmul(w, x, y, prec, rnd);
    }
    else if (arf_is_special(x) || arf_is_special(y) || arf_is_special(z))
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
            arf_mul(w, x, y, ARF_PREC_EXACT, ARF_RND_DOWN);
            return arf_add(w, z, w, prec, rnd);
        }
    }

    mp_size_t xn, yn, zn, tn, alloc;
    mp_srcptr xptr, yptr, zptr;
    mp_ptr tptr, tptr2;
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
    ARF_MUL_TMP_ALLOC(tptr2, alloc);
    tptr = tptr2;

    ARF_MPN_MUL(tptr, xptr, xn, yptr, yn);

    tn   -= (tptr[0] == 0);
    tptr += (tptr[0] == 0);

    if (shift >= 0)
        inexact = _arf_add_mpn(w, zptr, zn, ARF_SGNBIT(z), ARF_EXPREF(z),
            tptr, tn, tsgnbit, shift, prec, rnd);
    else
        inexact = _arf_add_mpn(w, tptr, tn, tsgnbit, texp,
            zptr, zn, ARF_SGNBIT(z), -shift, prec, rnd);

    ARF_MUL_TMP_FREE(tptr2, alloc)
    fmpz_clear(texp);

    return inexact;
}
