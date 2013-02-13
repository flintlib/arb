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
fmpr_get_si(const fmpr_t x, fmpr_rnd_t rnd)
{
    fmpz_t t;
    long v;
    fmpz_init(t);
    fmpr_get_fmpz(t, x, rnd);
    if (!fmpz_fits_si(t))
    {
        printf("fmpr_get_si: result does not fit in a signed long\n");
        abort();
    }
    v = fmpz_get_si(t);
    fmpz_clear(t);
    return v;
}

