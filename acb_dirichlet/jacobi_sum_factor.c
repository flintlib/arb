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
acb_dirichlet_jacobi_sum_factor(acb_t res,  const dirichlet_group_t G, const dirichlet_char_t chi1, const dirichlet_char_t chi2, slong prec)
{
    slong k;
    acb_t tmp;
    acb_init(tmp);
    acb_one(res);

    /* TODO: efficient subgroup */
    for (k = 0; k < G->num; k++)
    {
        nmod_t pe;
        ulong p, e, ap, bp;

        p = G->P[k].p;
        e = G->P[k].e;
        pe = G->P[k].pe;
        ap = chi1->x->n % pe.n;
        bp = chi2->x->n % pe.n;

        if (ap == 1 || bp == 1 || nmod_mul(ap, bp, pe) == 1)
        {
            slong r;
            ulong cond;

            cond = (ap == 1) ? chi2->conductor : chi1->conductor;
            r = jacobi_one_prime(p, e, pe.n, cond);

            /* chi(a,-1) if ap * bp = 1 */
            if (ap != 1 && bp != 1)
                r *= n_jacobi_unsigned(ap, p);

            acb_mul_si(res, res, r, prec);
        }
        else
        {
            dirichlet_group_t Gp;
            dirichlet_char_t chi1p, chi2p;

            dirichlet_group_init(Gp, pe.n);
            dirichlet_char_init(chi1p, Gp);
            dirichlet_char_init(chi2p, Gp);

            chi1p->x->n = ap;
            chi1p->x->log[0] = chi1->x->log[k];
            chi2p->x->n = ap;
            chi2p->x->log[0] = chi2->x->log[k];

            dirichlet_char_conrey(chi1p, Gp, NULL);
            dirichlet_char_conrey(chi2p, Gp, NULL);

            /* TODO: work out gauss relations for e > 1 */
            if (p <= 100 || e > 1)
                acb_dirichlet_jacobi_sum_naive(tmp, Gp, chi1p, chi2p, prec);
            else
                acb_dirichlet_jacobi_sum_gauss(tmp, Gp, chi1p, chi2p, prec);

            acb_mul(res, res, tmp, prec);

            dirichlet_char_clear(chi1p);
            dirichlet_char_clear(chi2p);
            dirichlet_group_clear(Gp);
        }
    }
    acb_clear(tmp);
}
