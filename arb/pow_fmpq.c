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

    Copyright (C) 2012, 2013 Fredrik Johansson

******************************************************************************/

#include "arb.h"

void
arb_pow_fmpq(arb_t y, const arb_t x, const fmpq_t a, long prec)
{
    if (fmpz_is_one(fmpq_denref(a)))
    {
        arb_pow_fmpz(y, x, fmpq_numref(a), prec);
    }
    else if (fmpz_cmp_ui(fmpq_denref(a), 36) <= 0)
    {
        arb_root(y, x, *fmpq_denref(a), prec);
        arb_pow_fmpz(y, y, fmpq_numref(a), prec);
    }
    else
    {
        long wp;

        wp = prec + 10;

        arb_log(y, x, wp);
        arb_mul_fmpz(y, y, fmpq_numref(a), wp);
        arb_div_fmpz(y, y, fmpq_denref(a), wp);
        arb_exp(y, y, prec);
    }
}

