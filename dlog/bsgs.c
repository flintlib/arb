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

    Copyright (C) 2016 Pascal Molin

******************************************************************************/

#include <stdlib.h>
#include "dlog.h"

ulong
dlog_bsgs(const dlog_bsgs_t t, ulong b)
{
    ulong i;
    apow_t c, * x;

    c.ak = b;
    for (i = 0; i < t->g; i++)
    {
        x = bsearch(&c, t->table, t->m, sizeof(apow_t),
            (int(*)(const void*,const void*))apow_cmp);
        if (x != NULL)
            return i * t->m + x->k;
        c.ak = nmod_mul(c.ak, t->am, t->mod);
    }
    flint_printf("Exception (n_discrete_log_bsgs).  discrete log not found.\n");
    flint_printf("   table size %wu, cosize %wu mod %wu. %wu not found (a^-m=%wu)\n",
            t->m, t->g, t->mod.n, b, t->am);
    abort();
}
