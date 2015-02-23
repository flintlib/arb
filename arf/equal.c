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

#include "arf.h"

int
arf_equal(const arf_t x, const arf_t y)
{
    mp_size_t n;

    if (x == y)
        return 1;

    if (ARF_XSIZE(x) != ARF_XSIZE(y))
        return 0;

    if (!fmpz_equal(ARF_EXPREF(x), ARF_EXPREF(y)))
        return 0;

    n = ARF_SIZE(x);

    if (n == 0)
        return 1;

    if (n == 1)
        return ARF_NOPTR_D(x)[0] == ARF_NOPTR_D(y)[0];

    if (n == 2)
        return (ARF_NOPTR_D(x)[0] == ARF_NOPTR_D(y)[0] &&
                ARF_NOPTR_D(x)[1] == ARF_NOPTR_D(y)[1]);

    return mpn_cmp(ARF_PTR_D(x), ARF_PTR_D(y), n) == 0;
}

int
arf_equal_si(const arf_t x, long y)
{
    arf_t t;
    arf_init_set_si(t, y);
    return arf_equal(x, t); /* no need to free */
}

