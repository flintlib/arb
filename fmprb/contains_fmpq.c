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

#include "fmprb.h"

int
fmprb_contains_fmpq(const fmprb_t x, const fmpq_t y)
{
    if (fmpz_is_one(fmpq_denref(y)))
    {
        return fmprb_contains_fmpz(x, fmpq_numref(y));
    }
    else
    {
        int ans;
        fmprb_t t;
        fmprb_init(t);
        fmprb_mul_fmpz(t, x, fmpq_denref(y), FMPR_PREC_EXACT);
        ans = fmprb_contains_fmpz(t, fmpq_numref(y));
        fmprb_clear(t);
        return ans;
    }
}
