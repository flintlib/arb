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
arb_pow_fmpq(arb_t y, const arb_t x, const fmpq_t a, slong prec)
{
    if (fmpz_is_one(fmpq_denref(a)))
    {
        arb_pow_fmpz(y, x, fmpq_numref(a), prec);
    }
    else
    {
        int use_exp;
        slong k = *fmpq_denref(a);

        if (k == 2 || k == 4)
            use_exp = 0;
        else if (k > 1 && k < 50)
            use_exp = prec < (1L << ((k / 8) + 8));
        else
            use_exp = 1;

        if (use_exp)
        {
            arb_log(y, x, prec + 10);
            arb_mul_fmpz(y, y, fmpq_numref(a), prec + 10);
            arb_div_fmpz(y, y, fmpq_denref(a), prec + 10);
            arb_exp(y, y, prec);
        }
        else
        {
            arb_root(y, x, k, prec);
            arb_pow_fmpz(y, y, fmpq_numref(a), prec);
        }
    }
}

