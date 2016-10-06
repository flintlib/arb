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
acb_dirichlet_dft_bluestein_init(acb_dirichlet_dft_bluestein_t t, slong n, slong prec)
{

    nmod_t n2;
    slong k, k2;
    acb_ptr z2n;
    int e = n_clog(2 * n - 1, 2);
    acb_dirichlet_dft_rad2_init(t->rad2, e, prec);
    z2n = _acb_vec_init(2 * n);
    acb_dirichlet_vec_nth_roots(z2n, 2 * n, prec);
    nmod_init(&n2, 2 * n);
    t->n = n;
    t->z = _acb_vec_init(n);
    for (k = 0, k2 = 0; k < n; k++)
    {
        acb_conj(t->z + k, z2n + k2);
        k2 = nmod_add(k2, 2 * k + 1, n2);
    }
    _acb_vec_clear(z2n, 2 * n);
}

void
acb_dirichlet_dft_bluestein_precomp(acb_ptr w, acb_srcptr v, const acb_dirichlet_dft_bluestein_t t, slong prec)
{
    slong n = t->n;
    acb_ptr vz, wz, z;
    z = t->z;
    /* TODO: allocate directly length 2^e and pad */
    flint_printf("\n\n====================\n\nv\n");
    acb_vec_printd_index(v, n, 10);
    vz = _acb_vec_init(n);
    acb_vec_kronecker_mul_conj(vz, z, v, n, prec);
    flint_printf("\nvz\n");
    acb_vec_printd_index(vz, n, 10);
    wz = _acb_vec_init(n);
    acb_dirichlet_dft_convol_rad2_precomp(wz, vz, z, n, t->rad2, prec);
    flint_printf("\nwz\n");
    acb_vec_printd_index(wz, n, 10);
    acb_vec_kronecker_mul_conj(w, z, wz, n, prec);
    flint_printf("\nw\n");
    acb_vec_printd_index(w, n, 10);
    _acb_vec_clear(wz, n);
    _acb_vec_clear(vz, n);
}

void
acb_dirichlet_dft_bluestein(acb_ptr w, acb_srcptr v, slong len, slong prec)
{
    acb_dirichlet_dft_bluestein_t t;
    acb_dirichlet_dft_bluestein_init(t, len, prec);
    acb_dirichlet_dft_bluestein_precomp(w, v, t, prec);
    acb_dirichlet_dft_bluestein_clear(t);
}
