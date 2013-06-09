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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "fmpr.h"

void
fmpr_divappr_abs_ubound(fmpr_t z, const fmpr_t x, const fmpr_t y, long prec)
{
    if (fmpr_is_special(x) || fmpr_is_special(y) || fmpz_is_pm1(fmpr_manref(y)))
    {
        fmpr_div(z, x, y, prec, FMPR_RND_UP);
        fmpr_abs(z, z);
    }
    else
    {
        fmpz_t t, u;
        long xbits, ybits, tbits, ubits, shift;

        xbits = fmpz_bits(fmpr_manref(x));
        ybits = fmpz_bits(fmpr_manref(y));

        fmpz_init(t);
        fmpz_init(u);

        ubits = FLINT_MIN(ybits, prec);
        tbits = prec + ubits + 1;

        /* upper bound for |x|, shifted */
        if (xbits <= tbits)
        {
            fmpz_mul_2exp(t, fmpr_manref(x), tbits - xbits);
            fmpz_abs(t, t);
        }
        else if (fmpz_sgn(fmpr_manref(x)) > 0)
        {
            fmpz_cdiv_q_2exp(t, fmpr_manref(x), xbits - tbits);
        }
        else
        {
            fmpz_fdiv_q_2exp(t, fmpr_manref(x), xbits - tbits);
            fmpz_neg(t, t);
        }

        /* lower bound for |y|, shifted */
        if (ybits <= ubits)
            fmpz_mul_2exp(u, fmpr_manref(y), ubits - ybits);
        else
            fmpz_tdiv_q_2exp(u, fmpr_manref(y), ybits - ubits);
        fmpz_abs(u, u);

        fmpz_cdiv_q(fmpr_manref(z), t, u);

        shift = (ubits - ybits) - (tbits - xbits);
        fmpz_sub(fmpr_expref(z), fmpr_expref(x), fmpr_expref(y));
        if (shift >= 0)
            fmpz_add_ui(fmpr_expref(z), fmpr_expref(z), shift);
        else
            fmpz_sub_ui(fmpr_expref(z), fmpr_expref(z), -shift);

        _fmpr_normalise(fmpr_manref(z), fmpr_expref(z), prec, FMPR_RND_UP);

        fmpz_clear(t);
        fmpz_clear(u);
    }
}

