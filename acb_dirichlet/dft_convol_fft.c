/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

/* assume np >= 2 * n - 1 */
void
acb_dirichlet_dft_convol_pad(acb_ptr fp, acb_ptr gp, acb_srcptr f, acb_srcptr g, slong n, slong np)
{
    slong k;
    for (k = 0; k < n; k++)
        acb_set(gp + k, g + k);
    for (k = 0; k < n; k++)
        acb_set(fp + k, f + k);
    for (k = 1; k < n; k++)
        acb_set(fp + np - k, f + n - k);
    for (k = n; k <= np - n; k++)
        acb_zero(fp + k);
}

void
acb_dirichlet_dft_inverse_cyc(acb_ptr w, acb_srcptr v, slong len, slong prec)
{
    /* divide before to keep v const */
    _acb_vec_scalar_div_ui(w, v, len, len, prec);
    acb_vec_swap_ri(w, len);
    acb_dirichlet_dft_cyc(w, w, len, prec);
    acb_vec_swap_ri(w, len);
}

void
acb_dirichlet_dft_convol_fft(acb_ptr w, acb_srcptr f, acb_srcptr g, slong len, slong prec)
{
    int e;
    slong k, np;
    acb_ptr fp, gp;
    acb_dirichlet_dft_rad2_t dft;
    e = n_clog(2 * len + 1, 2);
    acb_dirichlet_dft_rad2_init(dft, e, prec);
    np = dft->n;
    fp = _acb_vec_init(np);
    gp = _acb_vec_init(np);
    acb_dirichlet_dft_convol_pad(fp, gp, f, g, len, np);
    acb_dirichlet_dft_rad2_precomp(fp, dft, prec);
    acb_dirichlet_dft_rad2_precomp(gp, dft, prec);
    _acb_vec_kronecker_mul(gp, gp, fp, np, prec);
    acb_dirichlet_dft_inverse_rad2_precomp(gp, dft, prec);
    for (k = 0; k < len; k++)
        acb_set(w + k, gp + k);
}
