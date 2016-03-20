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

void
dlog_vec_set_not_found(ulong *v, ulong nv, nmod_t mod)
{
    n_factor_t fac;
    ulong i;

    n_factor_init(&fac);
    n_factor(&fac, mod.n, 1);
    for (i = 0; i < fac.num; i++)
    {
        ulong p, k;
        p = fac.p[i];
        for (k = p; k < nv; k += p)
            v[k] = NOT_FOUND;
    }
}
