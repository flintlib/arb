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

    Copyright (C) 2015 Fredrik Johansson
    Copyright (C) 2015 Arb authors

******************************************************************************/

#include "arb.h"

void
arb_fprint(FILE * file, const arb_t x)
{
    arf_fprint(file, arb_midref(x));
    flint_fprintf(file, " +/- ");
    mag_fprint(file, arb_radref(x));
}

void
arb_fprintd(FILE * file, const arb_t x, slong digits)
{
    arf_fprintd(file, arb_midref(x), FLINT_MAX(digits, 1));
    flint_fprintf(file, " +/- ");
    mag_fprintd(file, arb_radref(x), 5);
}

void
arb_fprintn(FILE * file, const arb_t x, slong digits, ulong flags)
{
    char * s = arb_get_str(x, digits, flags);
    flint_fprintf(file, "%s", s);
    flint_free(s);
}

