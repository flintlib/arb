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

    Copyright (C) 2015 Arb authors

******************************************************************************/

#include "fmpr.h"

void
fmpr_fprint(FILE * file, const fmpr_t x)
{
    if (fmpr_is_normal(x))
    {
        flint_fprintf(file, "(");
        fmpz_fprint(file, fmpr_manref(x));
        flint_fprintf(file, " * 2^");
        fmpz_fprint(file, fmpr_expref(x));
        flint_fprintf(file, ")");
    }
    else
    {
        if (fmpr_is_zero(x)) flint_fprintf(file, "(0)");
        else if (fmpr_is_pos_inf(x)) flint_fprintf(file, "(+inf)");
        else if (fmpr_is_neg_inf(x)) flint_fprintf(file, "(-inf)");
        else flint_fprintf(file, "(nan)");
    }
}
