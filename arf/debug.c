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

void
arf_debug(const arf_t x)
{
    mp_srcptr d;
    mp_size_t n;
    long i;

    printf("{exp=");
    fmpz_print(&x->exp);
    printf("; size=%lu; sgnbit=%ld; digits=[", ARF_SIZE(x), ARF_SGNBIT(x));

    ARF_GET_MPN_READONLY(d, n, x);

    for (i = 0; i < n; i++)
        printf(" %lu", d[i]);

    printf("]}");
}

