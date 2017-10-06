/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dft.h"

void
_acb_dft_bluestein_init(acb_dft_bluestein_t t, slong dv, slong n, slong prec)
{
    nmod_t n2;
    slong k, k2;
    acb_ptr z2n;
    int e = n_clog(2 * n - 1, 2);
#if DFT_VERB
    flint_printf("\n init bluestein e = %i", e);
#endif
    acb_dft_rad2_init(t->rad2, e, prec);

    /* compute z[k] = e(-k^2/2n) */
    /* TODO: check if this can be improved */
    z2n = _acb_vec_init(2 * n);
    _acb_vec_unit_roots(z2n, 2 * n, prec);
    nmod_init(&n2, 2 * n);
    t->n = n;
    t->dv = dv;
    t->z = _acb_vec_init(n);
    for (k = 0, k2 = 0; k < n; k++)
    {
        acb_set(t->z + k, z2n + k2);
        k2 = nmod_add(k2, 2 * k + 1, n2);
    }

    _acb_vec_clear(z2n, 2 * n);
}

void
acb_dft_bluestein_precomp(acb_ptr w, acb_srcptr v, const acb_dft_bluestein_t t, slong prec)
{
    slong k, n = t->n, np = t->rad2->n, dv = t->dv;
    acb_ptr fp, gp, z;
    z = t->z;

    fp = _acb_vec_init(np);
    _acb_vec_kronecker_mul_step(fp, z, v, dv, n, prec);

    gp = _acb_vec_init(np);
    acb_one(gp + 0);
    for (k = 1; k < n; k++)
    {
        acb_conj(gp + k, z + k);
        acb_set(gp + np - k, gp + k);
    }

    acb_dft_rad2_precomp_inplace(fp, t->rad2, prec);
    acb_dft_rad2_precomp_inplace(gp, t->rad2, prec);

    _acb_vec_kronecker_mul(gp, gp, fp, np, prec);

    acb_dft_inverse_rad2_precomp_inplace(gp, t->rad2, prec);

    _acb_vec_kronecker_mul(w, z, gp, n, prec);

    _acb_vec_clear(fp, n);
    _acb_vec_clear(gp, n);
}

void
acb_dft_bluestein(acb_ptr w, acb_srcptr v, slong len, slong prec)
{
    acb_dft_bluestein_t t;
    acb_dft_bluestein_init(t, len, prec);
    acb_dft_bluestein_precomp(w, v, t, prec);
    acb_dft_bluestein_clear(t);
}
