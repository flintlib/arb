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
    t->pe1 = pe / p;
    nmod_init(&t->pe, pe);
    t->inva = nmod_inv(a, t->pe);

    t->modp = flint_malloc(sizeof(dlog_precomp_struct));
    dlog_precomp_n_init(t->modp, a, p, p - 1, num);

    a1 = nmod_pow_ui(a, p - 1, t->pe);
    if (1 || e <= 2)
        dlog_1modpe_rec_init(t->modpe.rec, a1, p, e, t->pe);
    else
        dlog_1modpe_padic_init(t->modpe.padic, a1, p, e);

    return t->modp->cost + e;
}
