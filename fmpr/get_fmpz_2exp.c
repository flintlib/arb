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
fmpr_get_fmpz_2exp(fmpz_t man, fmpz_t exp, const fmpr_t x)
{
    if (fmpr_is_zero(x))
    {
        fmpz_zero(man);
        fmpz_zero(exp);
    }
    else
    {
        fmpz_set(man, fmpr_manref(x));
        fmpz_set(exp, fmpr_expref(x));
    }
}
