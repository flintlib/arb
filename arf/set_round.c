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
        /* XXX: fixme */
        if (y == x)
        {
            int inexact;
            arf_t t;
            arf_init(t);
            arf_set(t, x);
            inexact = arf_set_round(y, t, prec, rnd);
            arf_clear(t);
            return inexact;
        }
        else
        {
            mp_srcptr xptr;
            mp_size_t xn;
            long fix;
            int inexact;

            ARF_GET_MPN_READONLY(xptr, xn, x);
            inexact = _arf_set_round_mpn(y, &fix, xptr, xn,
                ARF_SGNBIT(x), prec, rnd);
            _fmpz_add_fast(ARF_EXPREF(y), ARF_EXPREF(x), fix);
            return inexact;
        }
    }
}

