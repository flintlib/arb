/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "arb.h"

void
_arb_rad_add_ufloat(arb_t y, const ufloat_t err)
{
    ufloat_t t, w;

    if (fmpz_is_zero(arb_radref(y)))
    {
        long yexp, eexp;

        yexp = *arb_expref(y);
        eexp = err->exp;

        if (yexp >= eexp)
        {
            fmpz_mul_2exp(arb_midref(y), arb_midref(y), yexp - eexp);
            fmpz_set_ui(arb_radref(y), err->man);
        }
        else
        {
            fmpz_tdiv_q_2exp(arb_midref(y), arb_midref(y), eexp - yexp);
            fmpz_set_ui(arb_radref(y), err->man + 1);
        }

        fmpz_set_si(arb_expref(y), eexp);
    }
    else
    {
        t->man = err->man;
        t->exp = err->exp - *arb_expref(y);

        ufloat_set_fmpz(w, arb_radref(y));
        ufloat_add(w, w, t);
        ufloat_get_fmpz(arb_radref(y), w);
    }
}
