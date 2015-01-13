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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "acb.h"

#ifndef __compar_fn_t
typedef int (*__compar_fn_t) (__const void *, __const void *);
#endif

int acb_cmp_pretty(const acb_t a, const acb_t b)
{
    arb_t t, u, v;
    int res;
    arb_init(t);
    arb_init(u);
    arb_init(v);
    arb_abs(u, acb_imagref(a));
    arb_abs(v, acb_imagref(b));
    arb_sub(t, u, v, MAG_BITS);
    res = 0;
    if (arb_contains_zero(t))
    {
        arb_sub(t, acb_realref(a), acb_realref(b), MAG_BITS);
        res = arb_is_positive(t) ? 1 : -1;
    }
    else
    {
        res = arb_is_positive(t) ? 1 : -1;
    }
    arb_clear(t);
    arb_clear(u);
    arb_clear(v);
    return res;
}

void _acb_vec_sort_pretty(acb_ptr vec, long len)
{
    qsort(vec, len, sizeof(acb_struct), (__compar_fn_t) acb_cmp_pretty);
}

