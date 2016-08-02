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

/* shallow copy */
static acb_ptr
vec_extract(acb_srcptr v, slong step, slong len)
{
    slong k;
    acb_ptr res;
    res = flint_malloc(len * sizeof(acb_struct));
    for (k = 0; k < len; k++)
    {
        res[k] = v[0];
        v += step;
    }
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
    slong num;
    slong n;
    slong * m;
    slong * M;
    acb_ptr z;
}
acb_dft_struct;

typedef acb_dft_struct acb_dft_t[1];

void
_acb_dft_fast(acb_ptr w, acb_srcptr v, slong dv, acb_srcptr z, slong dz, slong * cyc, slong * cocyc, slong num, slong prec)
{
    if (num == 1)
    {
        _acb_dft_base(w, v, dv, z, dz, cyc[0], prec);
    }
    else
    {
        slong i, j, m, M;
        acb_ptr t, wi;
        m = cyc[0];
        M = cocyc[0];
        t = _acb_vec_init(m * M);
        wi = w;
        for (i = 0; i < m; i++)
        {
            _acb_dft_fast(wi, v + i * dv, m * dv, z, dz * m, cyc + 1, cocyc + 1, num - 1, prec);
            if (i)
            { 
                for (j = 1; j < M; j++)
                    acb_mul(wi + j, wi + j, z + dz * i * j, prec);
            }
            wi += M;
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
acb_dirichlet_dft_fast(acb_ptr w, acb_srcptr v, slong len, slong prec)
{
    slong i, j, num;
    n_factor_t fac;
    slong * cyc, * cocyc;
    acb_ptr z;

    z = acb_roots_init(len, prec);

    n_factor_init(&fac);
    n_factor(&fac, len, 0);
    num = 0;
    for (i = 0; i < fac.num; i++)
        num += fac.exp[i];
    cyc = flint_malloc(num * sizeof(slong));
    cocyc = flint_malloc(num * sizeof(slong));

    num = 0;
    for (i = 0; i < fac.num; i++)
    {
        for (j = 0; j < fac.exp[i]; j++)
        {
            cyc[num] = fac.p[i];
            cocyc[num] = (len /= cyc[num]);
            num++;
        }
    }

    _acb_dft_fast(w, v, 1, z, 1, cyc, cocyc, num, prec);
}
