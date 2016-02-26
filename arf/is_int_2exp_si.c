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

    Copyright (C) 2016 Fredrik Johansson

******************************************************************************/

#include "arf.h"

int
arf_is_int_2exp_si(const arf_t x, slong e)
{
    if (arf_is_special(x))
        return arf_is_zero(x);
    else
    {
        fmpz_t t;
        int r;
        fmpz_init(t);
        arf_bot(t, x);
        r = fmpz_cmp_si(t, e) >= 0;
        fmpz_clear(t);
        return r;
    }
}

