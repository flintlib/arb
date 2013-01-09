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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "fmpr.h"

long
fmpr_set_fmpq(fmpr_t x, const fmpq_t y, long prec, fmpr_rnd_t rnd)
{
    if (fmpz_is_one(fmpq_denref(y)))
    {
        return fmpr_set_round_fmpz(x, fmpq_numref(y), prec, rnd);
    }
    else
    {
        long res;

        fmpr_t t, u;
        fmpr_init(t);
        fmpr_init(u);

        fmpr_set_fmpz(t, fmpq_numref(y));
        fmpr_set_fmpz(u, fmpq_denref(y));

        res = fmpr_div(x, t, u, prec, rnd);

        fmpr_clear(t);
        fmpr_clear(u);

        return res;
    }
}
