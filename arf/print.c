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

#include "arf.h"

void
arf_print(const arf_t x)
{
    if (arf_is_normal(x))
    {
        fmpz_t man, exp;

        fmpz_init(man);
        fmpz_init(exp);

        arf_get_fmpz_2exp(man, exp, x);

        printf("(");
        fmpz_print(man);
        printf(" * 2^");
        fmpz_print(exp);
        printf(")");

        fmpz_clear(man);
        fmpz_clear(exp);
    }
    else
    {
        if (arf_is_zero(x)) printf("(0)");
        else if (arf_is_pos_inf(x)) printf("(+inf)");
        else if (arf_is_neg_inf(x)) printf("(-inf)");
        else printf("(nan)");
    }
}

