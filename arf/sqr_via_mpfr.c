/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arf.h"

int
arf_sqr_via_mpfr(arf_ptr res, arf_srcptr x, slong prec, arf_rnd_t rnd)
{
    mp_size_t xn, resn, val;
    mp_srcptr xptr;
    mp_ptr tmp, resptr;
    mpfr_t xf, resf;
    int ret;
    ARF_MUL_TMP_DECL

    if (arf_is_special(x))
    {
        arf_sqr_special(res, x);
        return 0;
    }

    ARF_GET_MPN_READONLY(xptr, xn, x);

    prec = FLINT_MIN(xn * 2 * FLINT_BITS, prec);
    resn = (prec + FLINT_BITS - 1) / FLINT_BITS;

    ARF_MUL_TMP_ALLOC(tmp, resn)

    resf->_mpfr_d = tmp;
    resf->_mpfr_prec = prec;
    resf->_mpfr_sign = 1;
    resf->_mpfr_exp = 0;

    xf->_mpfr_d = (mp_ptr) xptr;
    xf->_mpfr_prec = xn * FLINT_BITS;
    xf->_mpfr_sign = 1; /* Sign does not matter since we are squaring */
    xf->_mpfr_exp = 0;

    ret = mpfr_sqr(resf, xf, arf_rnd_to_mpfr(rnd));

    ret = (ret != 0);

    _fmpz_2times_fast(ARF_EXPREF(res), ARF_EXPREF(x), resf->_mpfr_exp);

    val = 0;
    while (tmp[val] == 0)
        val++;

    ARF_GET_MPN_WRITE(resptr, resn - val, res);
    flint_mpn_copyi(resptr, tmp + val, resn - val);

    ARF_MUL_TMP_FREE(tmp, resn)

    return ret;
}
