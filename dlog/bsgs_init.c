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

    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2016 Pascal Molin

******************************************************************************/

#include <stdlib.h>
#include "dlog.h"

int
apow_cmp(const apow_t * x, const apow_t * y)
{
    return (x->ak < y->ak) ? -1 : (x->ak > y->ak);
}

/* set size of table m=sqrt(nk) to compute k logs in a group of size n */
ulong
dlog_bsgs_init(dlog_bsgs_t t, ulong a, ulong mod, ulong n, ulong m)
{
    ulong k, ak;
    if (m >= n) m = n + 1;
    t->table = (apow_t *)flint_malloc(m * sizeof(apow_t));

    nmod_init(&t->mod, mod);
    t->m = m;
    t->g = n / m + 1;

    for (k = 0, ak = 1; k < m; k++)
    {
        t->table[k].k = k;
        t->table[k].ak = ak;
        ak = nmod_mul(ak, a, t->mod);
    }

    t->am = nmod_inv(ak, t->mod);
    qsort(t->table, m, sizeof(apow_t), (int(*)(const void*,const void*))apow_cmp);
    return t->g;
}

void
dlog_bsgs_clear(dlog_bsgs_t t)
{
    flint_free(t->table);
}
