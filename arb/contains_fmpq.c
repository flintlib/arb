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

    Copyright (C) 2016 Fredrik Johansson

******************************************************************************/

#include "arb.h"

int
arb_contains_fmpq(const arb_t x, const fmpq_t y)
{
    if (fmpz_is_one(fmpq_denref(y)) || !arb_is_finite(x))
    {
        return arb_contains_fmpz(x, fmpq_numref(y));
    }
    else
    {
        arf_t t, xm, xr, ym;
        arf_struct tmp[3];
        int result;

        arf_init(t);
        arf_init(xm);
        arf_init(xr);
        arf_init(ym);

        /* To compare x with p/q, compare qx with p. */
        arf_mul_fmpz(xm, arb_midref(x), fmpq_denref(y), ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_set_mag(xr, arb_radref(x));
        arf_mul_fmpz(xr, xr, fmpq_denref(y), ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_set_fmpz(ym, fmpq_numref(y));

        /* y >= xm - xr  <=>  0 >= xm - xr - y */
        arf_init_set_shallow(tmp + 0, xm);
        arf_init_neg_shallow(tmp + 1,  xr);
        arf_init_neg_shallow(tmp + 2, ym);

        arf_sum(t, tmp, 3, 30, ARF_RND_DOWN);
        result = (arf_sgn(t) <= 0);

        if (result)
        {
            /* y <= xm + xr  <=>  0 <= xm + xr - y */
            arf_init_set_shallow(tmp + 1, xr);
            arf_sum(t, tmp, 3, 30, ARF_RND_DOWN);
            result = (arf_sgn(t) >= 0);
        }

        arf_clear(t);
        arf_clear(xm);
        arf_clear(xr);
        arf_clear(ym);

        return result;
    }
}

