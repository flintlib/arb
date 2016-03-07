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
#include "ulong_extras.h"

static int
apow_cmp(const apow_t * x, const apow_t * y)
{
    return (x->ak < y->ak) ? -1 : (x->ak > y->ak);
}

/* set size of table m=sqrt(nk) to compute k logs in a group of size n */
void
dlog_bsgs_init(dlog_bsgs_t t, ulong a, ulong n, ulong m)
{
    ulong k, ak;
    t->table = (apow_t *)flint_malloc(m * sizeof(apow_t));

    t->n = n;
    t->ninv = n_precompute_inverse(n);
    t->m = m;

    for (k = 0, ak = 1; k < m; k++)
    {
        t->table[k].k = k;
        t->table[k].ak = ak;
        ak = n_mulmod_precomp(ak, a, n, t->ninv);
    }

    t->am = n_invmod(ak, n);
    qsort(t->table, m, sizeof(apow_t), (int(*)(const void*,const void*))apow_cmp);
}

void
dlog_bsgs_clear(dlog_bsgs_t t)
{
    flint_free(t->table);
}

ulong
dlog_bsgs(const dlog_bsgs_t t, ulong b)
{
    ulong i;
    apow_t c, * x;

    c.k = 0;
    c.ak = b;
    for (i = 0; i < t->m; i++)
    {
        x = bsearch(&c, t->table, t->m, sizeof(apow_t),
            (int(*)(const void*,const void*))apow_cmp);
        if (x != NULL)
            return i * t->m + x->k;
        c.ak = n_mulmod_precomp(c.ak, t->am, t->n, t->ninv);
    }
    flint_printf("Exception (n_discrete_log_bsgs).  discrete log not found.\n");
    abort();
}
