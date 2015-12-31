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

#include "mag.h"
#include "arf.h"

void
mag_fprint(FILE * file, const mag_t x)
{
    flint_fprintf(file, "(");
    if (mag_is_zero(x))
        flint_fprintf(file, "0");
    else if (mag_is_inf(x))
        flint_fprintf(file, "inf");
    else
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_sub_ui(t, MAG_EXPREF(x), MAG_BITS);
        flint_fprintf(file, "%wu * 2^", MAG_MAN(x));
        fmpz_fprint(file, t);
        fmpz_clear(t);
    }
    flint_fprintf(file, ")");
}

void
mag_fprintd(FILE * file, const mag_t x, slong d)
{
    arf_t t;
    arf_init(t);
    arf_set_mag(t, x);
    arf_fprintd(file, t, d);
    arf_clear(t);
}
