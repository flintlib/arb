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
acb_dirichlet_gauss_sum_theta(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong prec)
{
    arb_t x;
    acb_t eps;

    arb_init(x);

    if (chi->conductor < G->q || (G->q == 300 && (chi->x->n == 71 || chi->x->n == 131))
            || (G->q == 600 && (chi->x->n == 11 || chi->x->n == 491)))
    {
        /* or could use l'Hopital rule */
        acb_dirichlet_gauss_sum_naive(res, G, chi, prec);
        return;
    }

    arb_one(x);
    acb_dirichlet_chi_theta_arb(res, G, chi, x, prec);
    acb_init(eps);
    acb_conj(eps, res);
    acb_div(eps, res, eps, prec);

    if (chi->parity)
    {
        arb_zero(acb_realref(res));
        arb_sqrt_ui(acb_imagref(res), G->q, prec);
    }
    else
    {
        arb_zero(acb_imagref(res));
        arb_sqrt_ui(acb_realref(res), G->q, prec);
    }

    acb_mul(res, res, eps, prec);

    arb_clear(x);
    acb_clear(eps);
}
