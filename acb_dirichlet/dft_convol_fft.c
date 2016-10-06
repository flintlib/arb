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

    if (np < 2 * n - 1)
    {
        flint_printf("dft_convol_pad: overlapping padding %ld < 2*%ld-1\n", np, n);
        abort();
    }

    for (k = 0; k < n; k++)
        acb_set(gp + k, g + k);
    for (; k < np; k++)
        acb_zero(gp + k);

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
acb_dirichlet_dft_convol_rad2_precomp(acb_ptr w, acb_srcptr f, acb_srcptr g, slong len, const acb_dirichlet_dft_rad2_t rad2, slong prec)
{
    slong np;
    acb_ptr fp, gp;
    np = rad2->n;

    flint_printf("\nf\n");
    acb_vec_printd_index(f, len, 10);
    flint_printf("\ng\n");
    acb_vec_printd_index(g, len, 10);

    fp = _acb_vec_init(np);
    gp = _acb_vec_init(np);
    acb_dirichlet_dft_convol_pad(fp, gp, f, g, len, np);

    flint_printf("\nF\n");
    acb_vec_printd_index(fp, np, 10);
    flint_printf("\nG\n");
    acb_vec_printd_index(gp, np, 10);

    acb_dirichlet_dft_rad2_precomp(fp, rad2, prec);

    flint_printf("\nDFT F\n");
    acb_vec_printd_index(fp, np, 10);

    acb_dirichlet_dft_rad2_precomp(gp, rad2, prec);

    flint_printf("\nDFT G\n");
    acb_vec_printd_index(gp, np, 10);

    _acb_vec_kronecker_mul(gp, gp, fp, np, prec);

    flint_printf("\n(DFT F)(DFT G)=DFT(F*G)\n");
    acb_vec_printd_index(gp, np, 10);

    acb_dirichlet_dft_inverse_rad2_precomp(gp, rad2, prec);

    flint_printf("\nF*G\n");
    acb_vec_printd_index(gp, np, 10);

    _acb_vec_set(w, gp, len);
    _acb_vec_clear(fp, np);
    _acb_vec_clear(gp, np);
}

void
acb_dirichlet_dft_convol_rad2(acb_ptr w, acb_srcptr f, acb_srcptr g, slong len, slong prec)
{
    int e;
    acb_dirichlet_dft_rad2_t dft;
    e = n_clog(2 * len - 1, 2);
    acb_dirichlet_dft_rad2_init(dft, e, prec);
    acb_dirichlet_dft_convol_rad2_precomp(w, f, g, len, dft, prec);
    acb_dirichlet_dft_rad2_clear(dft);
}
