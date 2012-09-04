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

void
fmpr_get_fmpq(fmpq_t y, const fmpr_t x)
{
    if (fmpr_is_zero(x))
    {
        fmpq_zero(y);
    }
    else if (fmpr_is_special(x) || COEFF_IS_MPZ(*fmpr_expref(x)))
    {
        printf("exception: fmpr_get_fmpq: cannot convert to rational\n");
        abort();
    }
    else
    {
        long exp = *fmpr_expref(x);

        fmpz_set_ui(fmpq_denref(y), 1UL);

        if (exp >= 0)
        {
            fmpz_mul_2exp(fmpq_numref(y), fmpr_manref(x), exp);
        }
        else
        {
            fmpz_set(fmpq_numref(y), fmpr_manref(x));
            fmpz_mul_2exp(fmpq_denref(y), fmpq_denref(y), -exp);
        }
    }
}
