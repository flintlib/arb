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
dlog_crt_init(dlog_crt_t t, ulong a, ulong mod, ulong n, ulong num)
{
    int k;
    n_factor_t fac;
    ulong * M, * u;
    ulong cost = 0;

    n_factor_init(&fac);
    n_factor(&fac, n, 1);

    t->num = fac.num;
    nmod_init(&t->mod,mod);
    nmod_init(&t->n, n);

    M = t->expo = flint_malloc(t->num * sizeof(ulong));
    u = t->crt_coeffs = flint_malloc(t->num * sizeof(ulong));
    t->pre = flint_malloc(t->num * sizeof(dlog_precomp_t));
    for (k = 0; k < t->num; k++)
        t->pre[k] = flint_malloc(sizeof(dlog_precomp_struct));

    for (k = 0; k < t->num; k++)
    {
        ulong p, e, mk;
        p = fac.p[k];
        e = fac.exp[k];
        mk = n_pow(p, e);
        M[k] = n / mk;
        u[k] = nmod_mul(M[k], n_invmod(M[k] % mk, mk), t->n);
        /* depends on the power */
#if 0
        flint_printf("[sub-crt -- init for size %wu mod %wu]\n", mk, mod);
#endif
        dlog_precomp_pe_init(t->pre[k], nmod_pow_ui(a, M[k], t->mod), mod, p, e, mk, num);
        cost += t->pre[k]->cost;
    }
#if 0
    if (cost > 500)
    flint_printf("[crt init for size %wu mod %wu -> cost %wu]\n", n,mod,cost);
#endif
    return cost;
}

void
dlog_crt_clear(dlog_crt_t t)
{
    int k;
    flint_free(t->expo);
    flint_free(t->crt_coeffs);
    for (k = 0; k < t->num; k++)
    {
        dlog_precomp_clear(t->pre[k]);
        flint_free(t->pre[k]);
    }
    flint_free(t->pre);
}
