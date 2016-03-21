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

/* vector of log(k,a)*loga % order in Z/modZ */ 
void
dlog_vec_loop(ulong * v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order)
{
    ulong x, vx;
    dlog_vec_fill(v, nv, DLOG_NOT_FOUND);
    x = 1; vx = 0;
    do
    {
        if (x < nv)
            v[x] = vx;
        x = nmod_mul(x, a, mod);
        vx = nmod_add(vx, va, order);
    }
    while (x != 1);
    for(x = mod.n + 1; x < nv; x++)
        v[x] = v[x - mod.n];
}
