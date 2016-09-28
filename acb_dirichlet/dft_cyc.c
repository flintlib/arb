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
_acb_dirichlet_dft_cyc_init_z_fac(acb_dirichlet_dft_cyc_t t, slong len, n_factor_t fac, slong dv, acb_ptr z, slong dz, slong prec)
{
    slong i, j, num;
    t->n = len;
    num = 0;
    for (i = 0; i < fac.num; i++)
        num += fac.exp[i];
    t->num = num;
    t->cyc = flint_malloc(num * sizeof(acb_dirichlet_dft_step_struct));

    /* FIXME: if z is allocated here, it must be cleared somewhere */
    if (z == NULL)
    {
        z = _acb_vec_init(t->n);
        acb_dirichlet_vec_nth_roots(z, t->n, prec);
        dz = 1;
    }
    t->z = z;

    num = 0;
    for (i = 0; i < fac.num; i++)
    {
        for (j = 0; j < fac.exp[i]; j++)
        {
            slong m, M;
            m = fac.p[i];
            M = (len /= m);
            t->cyc[num].m = m;
            t->cyc[num].M = M;
            t->cyc[num].dv = dv;
            t->cyc[num].z = z;
            t->cyc[num].dz = dz;
            _acb_dirichlet_dft_precomp_init(t->cyc[num].pre, M, z, dz * M, m, prec);
            dv *= m;
            dz *= m;
            num++;
        }
    }
 }

void
_acb_dirichlet_dft_cyc_init(acb_dirichlet_dft_cyc_t t, slong dv, slong len, slong prec)
{
    n_factor_t fac;
    n_factor_init(&fac);
    n_factor(&fac, len, 0);
    _acb_dirichlet_dft_cyc_init_z_fac(t, len, fac, dv, NULL, 0, prec);
}

void
acb_dirichlet_dft_cyc_clear(acb_dirichlet_dft_cyc_t t)
{
    slong i;
    for (i = 0; i < t->num; i++)
        acb_dirichlet_dft_precomp_clear(t->cyc[i].pre);
    flint_free(t->cyc);
    _acb_vec_clear(t->z, t->n);
}


#if 0
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
_acb_dirichlet_dft_base(acb_ptr w, acb_srcptr v, slong dv, acb_srcptr z, slong dz, slong n, slong prec)
{
    acb_ptr v1, z1;
    v1 = vec_extract(v, dv, n);
    z1 = vec_extract(z, dz, n);
    _acb_dirichlet_dft_pol(w, v1, z1, n, prec);
    flint_free(v1);
    flint_free(z1);
}
#endif

void
acb_dirichlet_dft_cyc(acb_ptr w, acb_srcptr v, slong len, slong prec)
{
    acb_dirichlet_dft_cyc_t t;
    acb_dirichlet_dft_cyc_init(t, len, prec);
    acb_dirichlet_dft_cyc_precomp(w, v, t, prec);
    acb_dirichlet_dft_cyc_clear(t);
}
