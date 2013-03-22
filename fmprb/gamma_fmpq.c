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

#include "fmprb.h"
#include "gamma.h"

void
fmprb_gamma_fmpq(fmprb_t y, const fmpq_t x, long prec)
{
    fmpz p, q;

    p = *fmpq_numref(x);
    q = *fmpq_denref(x);

    if ((q == 1 || q == 2 || q == 3 || q == 4 || q == 6) && !COEFF_IS_MPZ(p))
    {
        if (q == 1)
        {
            if (p <= 0)
            {
                fmpr_nan(fmprb_midref(y));
                fmpr_pos_inf(fmprb_radref(y));
                return;
            }

            if (p < 1200 || 1.44265 * (p*log(p) - p) < 15.0 * prec)
            {
                fmpz_t t;
                fmpz_init(t);
                fmpz_fac_ui(t, p - 1);
                fmprb_set_round_fmpz(y, t, prec);
                fmpz_clear(t);
                return;
            }
        }

        p = FLINT_ABS(p);

        if (p < q * 500.0 || p < q * (500.0 + 0.1 * prec * sqrt(prec)))
        {
            gamma_fmpq_outward(y, x, prec);
            return;
        }
    }

    gamma_fmpq_stirling(y, x, prec);
}

