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

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include "arf.h"

TLS_PREFIX mp_ptr __arf_mul_tmp = NULL;
TLS_PREFIX long __arf_mul_alloc = 0;

void _arf_mul_tmp_cleanup(void)
{
    flint_free(__arf_mul_tmp);
    __arf_mul_tmp = NULL;
    __arf_mul_alloc = 0;
}

int
arf_mul_rnd_any(arf_ptr z, arf_srcptr x, arf_srcptr y,
        long prec, arf_rnd_t rnd)
{
    mp_size_t xn, yn;
    long fix;
    int sgnbit, inexact;

    xn = ARF_XSIZE(x);
    yn = ARF_XSIZE(y);
    sgnbit = (xn ^ yn) & 1;
    xn >>= 1;
    yn >>= 1;

    if (yn > xn)
    {
        arf_srcptr __t; mp_size_t __u;
        __t = x; x = y; y = __t;
        __u = xn; xn = yn; yn = __u;
    }

    /* Either operand is a special value. */
    if (yn == 0)
    {
        arf_mul_special(z, x, y);
        return 0;
    }
    else
    {
        mp_size_t zn, alloc;
        mp_srcptr xptr, yptr;
        mp_ptr tmp;
        ARF_MUL_TMP_DECL

        ARF_GET_MPN_READONLY(xptr, xn, x);
        ARF_GET_MPN_READONLY(yptr, yn, y);

        alloc = zn = xn + yn;
        ARF_MUL_TMP_ALLOC(tmp, alloc)

        if (yn == 1)
        {
            tmp[zn - 1] = mpn_mul_1(tmp, xptr, xn, yptr[0]);
        }
        else if (xn == yn)
        {
            if (xptr == yptr)
                mpn_sqr(tmp, xptr, xn);
            else
                mpn_mul_n(tmp, xptr, yptr, yn);
        }
        else 
        {
            mpn_mul(tmp, xptr, xn, yptr, yn);
        }

        inexact = _arf_set_round_mpn(z, &fix, tmp, zn, sgnbit, prec, rnd);
        _fmpz_add2_fast(ARF_EXPREF(z), ARF_EXPREF(x), ARF_EXPREF(y), fix);
        ARF_MUL_TMP_FREE(tmp, alloc)

        return inexact;
    }
}

