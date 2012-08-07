/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "ufloat.h"

void
ufloat_log1p(ufloat_t z, const ufloat_t x)
{
    ufloat_t t;

    if (x->exp >= -1L)
    {
        ufloat_one(t);
        ufloat_add(t, t, x);
        ufloat_log(z, t);
    }
    else
    {
        /* log(1+x) <= x, todo: improve accuracy */
        ufloat_set(z, x);
    }
}
