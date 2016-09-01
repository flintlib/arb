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

#include "acb_dirichlet.h"

/* shallow extract */
static acb_ptr
vec_extract(acb_srcptr v, slong step, slong len)
{
    slong k, l;
    acb_ptr res;
    res = flint_malloc(len * sizeof(acb_struct));
    for (k = 0, l = 0; k < len; k++, l+=step)
        res[k] = v[l];
    return res;
}
void
_acb_dft_base(acb_ptr w, acb_srcptr v, slong dv, acb_srcptr z, slong dz, slong n, slong prec)
{
    acb_ptr v1, z1;
    v1 = vec_extract(v, dv, n);
    z1 = vec_extract(z, dz, n);
    _acb_dirichlet_dft_pol(w, v1, z1, n, prec);
    flint_free(v1);
    flint_free(z1);
}

typedef struct
{
    slong m;
    slong M;
    slong dv;
    acb_srcptr z;
    slong dz;
}
dft_cyc_step;

typedef struct
{
    slong n;
    acb_ptr z;
    slong num;
    dft_cyc_step * cyc;
}
acb_dft_cyc_struct;

typedef acb_dft_cyc_struct acb_dft_cyc_t[1];

typedef struct
{
    slong m;
    slong M;
    slong dv;
    acb_dft_cyc_t pre;
}
dft_prod_step;

typedef struct
{
    slong n;
    slong num;
    dft_prod_step * cyc;
}
acb_dft_prod_struct;

typedef acb_dft_prod_struct acb_dft_prod_t[1];

void
_acb_dirichlet_dft_cyc_init(acb_dft_cyc_t t, slong len, slong dv, slong prec)
{
    slong i, j, num, dz;
    n_factor_t fac;
    acb_ptr z;

    t->n = len;

    n_factor_init(&fac);
    n_factor(&fac, len, 0);
    num = 0;
    for (i = 0; i < fac.num; i++)
        num += fac.exp[i];
    t->num = num;
    t->cyc = flint_malloc(num * sizeof(dft_cyc_step));
    t->z = z = _acb_vec_init(t->n);
    acb_dirichlet_vec_nth_roots(z, t->n, prec);

    num = 0;
    dz = 1;
    for (i = 0; i < fac.num; i++)
    {
        for (j = 0; j < fac.exp[i]; j++)
        {
            slong m = fac.p[i];
            t->cyc[num].m = m;
            t->cyc[num].M = (len /= m);
            t->cyc[num].dv = dv;
            t->cyc[num].z = z;
            t->cyc[num].dz = dz;
            dv *= m;
            dz *= m;
            num++;
        }
    }
}

static void
acb_dirichlet_dft_cyc_init(acb_dft_cyc_t t, slong len, slong prec)
{
    _acb_dirichlet_dft_cyc_init(t, len, 1, prec);
}

void
acb_dirichlet_dft_cyc_clear(acb_dft_cyc_t t)
{
    _acb_vec_clear(t->z, t->n);
    flint_free(t->cyc);
}

void
_acb_dft_cyc(acb_ptr w, acb_srcptr v, dft_cyc_step * cyc, slong num, slong prec)
{
    dft_cyc_step c;
    if (num == 0)
        abort(); /* or just copy v to w */
    c = cyc[0];
    if (num == 1)
    {
        _acb_dft_base(w, v, c.dv, c.z, c.dz, c.m, prec);
    }
    else
    {
        slong i, j;
        slong m = c.m, M = c.M, dv = c.dv, dz = c.dz;
        acb_srcptr z = c.z, vi;
        acb_ptr t, wi;

        t = _acb_vec_init(m * M);
        wi = w; vi = v;
        for (i = 0; i < m; i++)
        {
            _acb_dft_cyc(wi, vi, cyc + 1, num - 1, prec);
            if (i)
            {
                for (j = 1; j < M; j++)
                    acb_mul(wi + j, wi + j, z + dz * i * j, prec);
            }
            wi += M;
            vi += dv;
        }
        /* after first pass */
        for (j = 0; j < M; j++)
            _acb_dft_base(t + m * j, w + j, M, z, dz * M, m, prec);
        /* reorder */
        for (i = 0; i < m; i++)
            for (j = 0; j < M; j++)
                w[j + M * i] = t[i + m * j];
        _acb_vec_clear(t, m * M);
    }
}

void
acb_dirichlet_dft_cyc_precomp(acb_ptr w, acb_srcptr v, acb_dft_cyc_t t, slong prec)
{
    _acb_dft_cyc(w, v, t->cyc, t->num, prec);
}

void
acb_dirichlet_dft_fast(acb_ptr w, acb_srcptr v, slong len, slong prec)
{
    acb_dft_cyc_t t;
    acb_dirichlet_dft_cyc_init(t, len, prec);
    acb_dirichlet_dft_cyc_precomp(w, v, t, prec);
    acb_dirichlet_dft_cyc_clear(t);
}

void
_acb_dft_prod(acb_ptr w, acb_srcptr v, dft_prod_step * cyc, slong num, slong prec)
{
    dft_prod_step c;
    if (num == 0)
    {
        flint_printf("error: reached num = 0 in dft_prod\n");
        abort(); /* or just copy v to w */
    }
    c = cyc[0];
    if (num == 1)
    {
        acb_dirichlet_dft_cyc_precomp(w, v, c.pre, prec);
    }
    else
    {
        slong i, j;
        slong m = c.m, M = c.M, dv = c.dv;
        acb_srcptr vi;
        acb_ptr t, wi;

        t = _acb_vec_init(m * M);
        wi = w; vi = v;
        for (i = 0; i < m; i++)
        {
            _acb_dft_prod(wi, vi, cyc + 1, num - 1, prec);
            wi += M;
            vi += dv; /* here = M */
        }
        /* after first pass */
        for (j = 0; j < M; j++)
            acb_dirichlet_dft_cyc_precomp(t + m * j, w + j, c.pre, prec);
        /* reorder */
        for (i = 0; i < m; i++)
            for (j = 0; j < M; j++)
                w[j + M * i] = t[i + m * j];
        _acb_vec_clear(t, m * M);
    }
}

void
acb_dirichlet_dft_prod_precomp(acb_ptr w, acb_srcptr v, acb_dft_prod_t t, slong prec)
{
    if (t->num == 0)
    {
        acb_set(w + 0, v + 0);
    }
    else
    {
        _acb_dft_prod(w, v, t->cyc, t->num, prec);
    }
}

void
acb_dirichlet_dft_prod_init(acb_dft_prod_t t, slong * cyc, slong num, slong prec)
{
    slong i, len, dv;

    t->num = num;
    t->cyc = flint_malloc(num * sizeof(dft_prod_step));

    len = 1; dv = 1;
    for (i = 0; i < num; i++)
        len *= cyc[i];

    for (i = 0; i < num; i++)
    {
        slong m = cyc[i];
        len /= m;
        t->cyc[i].m = m;
        t->cyc[i].M = len;
        t->cyc[i].dv = len;
        _acb_dirichlet_dft_cyc_init(t->cyc[i].pre, m, len, prec);
        dv *= m;
    }
}

void
acb_dirichlet_dft_prod_clear(acb_dft_prod_t t)
{
    slong i;
    for (i = 0; i < t->num; i++)
        acb_dirichlet_dft_cyc_clear(t->cyc[i].pre);
    flint_free(t->cyc);
}

void
acb_dirichlet_dft_prod(acb_ptr w, acb_srcptr v, slong * cyc, slong num, slong prec)
{
    acb_dft_prod_t t;
    acb_dirichlet_dft_prod_init(t, cyc, num, prec);
    acb_dirichlet_dft_prod_precomp(w, v, t, prec);
    acb_dirichlet_dft_prod_clear(t);
}
