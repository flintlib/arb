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

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include "fmpz_extras.h"

long
_fmpz_sub_small_large(const fmpz_t x, const fmpz_t y)
{
    fmpz_t t;
    fmpz_init(t);
    fmpz_sub(t, x, y);
    if (!COEFF_IS_MPZ(*t))
    {
        /* no need to free t */
        return *t;
    }
    else
    {
        int sign = fmpz_sgn(t);
        fmpz_clear(t);
        return (sign > 0) ? LONG_MAX : -LONG_MAX;
    }
}
