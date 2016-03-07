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
dlog_crt_init(dlog_crt_t t, ulong a, ulong mod, ulong n, ulong num)
{
    int k;
    n_factor_t fac;
    ulong * M, * u;

    n_factor_init(&fac);
    n_factor(&fac, n, 1);

    t->num = fac.num;
    t->mod = mod;
    t->n = n;

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
        u[k] = M[k] * n_invmod(M[k] % mk, mk) % n;
        /* depends on the power */
        dlog_precomp_pe_init(t->pre[k], n_powmod(a, M[k], mod), mod, p, e, mk, num);
    }
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

ulong
dlog_crt(const dlog_crt_t t, ulong b)
{
    int k;
    ulong r = 0;
    for (k = 0; k < t->num; k++)
    {
        ulong bk, rk;
        bk = n_powmod(b, t->expo[k], t->mod);
        rk = dlog_precomp(t->pre[k], bk);
#if 0
        flint_printf("##[crt-%d]: log(%wu)=log(%wu^%wu) = %wu [size %wu mod %wu]\n",
                k, bk, b, t->expo[k], rk, t->n/t->expo[k], t->mod);
#endif
        r = (r + t->crt_coeffs[k] * rk) % t->n;
    }
    return r;
}
