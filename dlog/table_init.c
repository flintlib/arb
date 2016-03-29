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

    Copyright (C) 2016 Pascal Molin

******************************************************************************/

#include "dlog.h"

/* assume mod is small so no overflow */
ulong
dlog_table_init(dlog_table_t t, ulong a, ulong mod)
{
    int k;
    ulong ak;
    t->mod = mod;
    t->table = flint_malloc(mod * sizeof(ulong));
    ak = 1; k = 0;

    /* warning: do not check a is invertible modulo mod */
    do
    {
        t->table[ak] = k++;
        ak = (ak * a) % mod;
    }

    while (ak != 1);
    return 1;
}
