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

int
arf_set_round(arf_t y, const arf_t x, long prec, arf_rnd_t rnd)
{
    if (arf_is_special(x))
    {
        arf_set(y, x);
        return 0;
    }
    else
    {
        int inexact;
        long fix;
        mp_size_t xn;
        mp_srcptr xptr;

        if (y == x)
        {
            mp_ptr xtmp;
            TMP_INIT;

            ARF_GET_MPN_READONLY(xptr, xn, x);

            /* exact */
            if (xn * FLINT_BITS <= prec)
                return 0;

            if ((xn - 1) * FLINT_BITS < prec)
            {
                /* exact */
                if ((xptr[0] << (prec - (xn-1) *  FLINT_BITS)) == 0)
                    return 0;
            }

            /* inexact */
            TMP_START;
            xtmp = TMP_ALLOC(xn * sizeof(mp_limb_t));
            flint_mpn_copyi(xtmp, xptr, xn);
            inexact = _arf_set_round_mpn(y, &fix, xtmp, xn, ARF_SGNBIT(x), prec, rnd);
            _fmpz_add_fast(ARF_EXPREF(y), ARF_EXPREF(x), fix);

            TMP_END;
            return inexact;
        }
        else
        {
            ARF_GET_MPN_READONLY(xptr, xn, x);
            inexact = _arf_set_round_mpn(y, &fix, xptr, xn,
                ARF_SGNBIT(x), prec, rnd);
            _fmpz_add_fast(ARF_EXPREF(y), ARF_EXPREF(x), fix);
            return inexact;
        }
    }
}

