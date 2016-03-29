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

ulong
dlog_modpe_init(dlog_modpe_t t, ulong a, ulong p, ulong e, ulong pe, ulong num)
{
    ulong a1;

    t->p = p;
    t->e = e;
    nmod_init(&t->pe, pe);
    t->inva = nmod_inv(a, t->pe);

    if (p == 2)
    {
        t->modp = NULL;
        t->pe1 = (e <= 2) ? 2 : pe / 4;
        t->modpe->inv1p = t->inva;
        t->modpe->invloga1 = 1;
        return e - 2;
    }
    else
    {
        t->modp = flint_malloc(sizeof(dlog_precomp_struct));
        t->pe1 = pe / p;
        dlog_precomp_n_init(t->modp, a, p, p - 1, num);

        a1 = nmod_pow_ui(a, p - 1, t->pe);
        dlog_1modpe_init(t->modpe, a1, p, e, t->pe);

        return t->modp->cost + e;
    }
}
