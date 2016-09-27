/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

void
crt_init(crt_t c, ulong n)
{
    int k;
    n_factor_t fac;

    n_factor_init(&fac);
    n_factor(&fac, n, 1);

    nmod_init(&c->n, n);

    c->num = fac.num;
    for (k = 0; k < fac.num; k++)
    {
        c->m[k] = n_pow(fac.p[k], fac.exp[k]);
        c->M[k] = n / c->m[k];
        c->vM[k] = c->M[k] * n_invmod(c->M[k] % c->m[k], c->m[k]);
        /*
        flint_printf("m[%ld]=%ld, M[%ld]=%wu, vM[%ld]=%wu\n",
                k, c->m[k], k, c->M[k], k, c->vM[k]);
                */
    }
}

#if 0
/* lexicographic index of crt elt j */
static ulong
index_crt(crt_t c, ulong j)
{
    int k;
    ulong res = 0;
    for (k = 0; k < c->num; k ++)
        res = res * c->m[k] + (j % c->m[k]);
    return res;
}

/* crt elt of lexicographic index i */
static ulong
crt_index(crt_t c, ulong i)
{
    int k;
    ulong j, res = 0;
    for (k = 0; k < c->num; k ++)
    {
        j = i % c->m[k];
        i = i / c->m[k];
        res = nmod_add(res,  j * c->vM[k], c->n);
    }
    return res;
}
/* for all elements can do fast conrey-like loop just adding vM[k] */
static acb_ptr
reorder_to_crt(acb_srcptr v, crt_t c, ulong len)
{
    ulong k;
    acb_ptr res;
    res = flint_malloc(len * sizeof(acb_struct));
    for (k = 0; k < len; k++)
        res[index_crt(c, k)] = v[k];
    return res;
}
static acb_ptr
reorder_from_crt(acb_srcptr v, crt_t c, ulong len)
{
    ulong k;
    acb_ptr res;
    res = flint_malloc(len * sizeof(acb_struct));
    for (k = 0; k < len; k++)
        res[k] = v[index_crt(c, k)];
    return res;
}
#endif

void
crt_decomp(acb_ptr y, acb_srcptr x, const crt_t c, ulong len)
{
    int j, e[CRT_MAX];
    ulong k, l;

    for (j = 0; j < c->num; j++)
        e[j] = 0;

    l = 0;
    for(k = 0; k < len; k++)
    {
        /*flint_printf("set y[%wu] = x[%wu]\n", k, l);*/
        acb_set(y + k, x + l);
        for (j = c->num - 1; j >= 0; e[j] = 0, j--)
        {
            e[j]++; l = nmod_add(l, c->vM[j], c->n);
            if (e[j] < c->m[j])
                break;
        }
    }
}

void
crt_recomp(acb_ptr y, acb_srcptr x, const crt_t c, ulong len)
{
    int j, e[CRT_MAX];
    ulong k, l;

    for (j = 0; j < c->num; j++)
        e[j] = 0;

    l = 0;
    for(k = 0; k < len; k++)
    {
        acb_set(y + l, x + k);
        for (j = c->num - 1; j >= 0; e[j] = 0, j--)
        {
            e[j]++; l = nmod_add(l, c->M[j], c->n);
            if (e[j] < c->m[j])
                break;
        }
    }
}

void
acb_dirichlet_dft_crt(acb_ptr w, acb_srcptr v, slong len, slong prec)
{
    crt_t c;
    acb_ptr t;
    t = _acb_vec_init(len);
    crt_init(c, len);
    crt_decomp(w, v, c, len);
    acb_dirichlet_dft_prod(t, w, c->m, c->num, prec);
    crt_recomp(w, t, c, len);
    _acb_vec_clear(t, len);
}
