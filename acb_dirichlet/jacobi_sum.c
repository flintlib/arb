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

/* J_N(1,a) = sum on x = 1 mod some p | q */
static ulong
jacobi_one_prime(ulong p, ulong e, ulong pe, ulong cond)
{
    if (e > 1 && cond % (p*p) == 0)
    {
        return 0;
    }
    else
    {
        slong r = (e > 1) ? pe / p : 1;
        if (cond % p)
            return r * (p - 2);
        else
            return -r;
    }
}

static ulong
jacobi_one(const acb_dirichlet_group_t G, ulong cond)
{
    slong k, r = 1;

    for (k = 0; k < G->num; k++)
        r *= jacobi_one_prime(G->P[k].p, G->P[k].e,
                G->P[k].pe.n, cond);
    return r;
}

/* should use only for prime power modulus */
static void
acb_dirichlet_jacobi_sum_gauss(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi1, const acb_dirichlet_char_t chi2, slong prec)
{
    /* J_q(a,b)G_q(ab) = G_q(a)G_q(b) */
    acb_t tmp;
    acb_dirichlet_char_t chi12;

    acb_dirichlet_char_init(chi12, G);
    acb_dirichlet_char_mul(chi12, G, chi1, chi2);

    acb_init(tmp);

    acb_dirichlet_gauss_sum(res, G, chi1, prec);
    if (chi2->x->n == chi1->x->n)
        acb_set(tmp, res);
    else
        acb_dirichlet_gauss_sum(tmp, G, chi2, prec);
    acb_mul(res, res, tmp, prec);
    acb_dirichlet_gauss_sum(tmp, G, chi12, prec);
    acb_div(res, res, tmp, prec);

    acb_dirichlet_char_clear(chi12);
    acb_clear(tmp);
}

static void
acb_dirichlet_jacobi_sum_primes(acb_t res,  const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi1, const acb_dirichlet_char_t chi2, slong prec)
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
            acb_dirichlet_group_t Gp;
            acb_dirichlet_char_t chi1p, chi2p;

            acb_dirichlet_group_init(Gp, pe.n);
            acb_dirichlet_char_init(chi1p, Gp);
            acb_dirichlet_char_init(chi2p, Gp);

            chi1p->x->n = ap;
            chi1p->x->log[0] = chi1->x->log[k];
            chi2p->x->n = ap;
            chi2p->x->log[0] = chi2->x->log[k];

            acb_dirichlet_char_conrey(chi1p, Gp, NULL);
            acb_dirichlet_char_conrey(chi2p, Gp, NULL);

            /* TODO: work out gauss relations for e > 1 */
            if (p <= 100 || e > 1)
                acb_dirichlet_jacobi_sum_naive(tmp, Gp, chi1p, chi2p, prec);
            else
                acb_dirichlet_jacobi_sum_gauss(tmp, Gp, chi1p, chi2p, prec);

            acb_mul(res, res, tmp, prec);

            acb_dirichlet_char_clear(chi1p);
            acb_dirichlet_char_clear(chi2p);
            acb_dirichlet_group_clear(Gp);
        }
    }
    acb_clear(tmp);
}

void
acb_dirichlet_jacobi_sum(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi1, const acb_dirichlet_char_t chi2, slong prec)
{

    if (G->q_even > 1)
    {
        acb_zero(res);
    }
    else if (chi1->x->n == 1 || chi2->x->n == 1)
    {
        ulong cond = (chi1->x->n == 1) ? chi2->conductor : chi1->conductor;
        acb_set_si(res, jacobi_one(G, cond));
    }
    else if (nmod_mul(chi1->x->n, chi2->x->n, G->mod) == 1)
    {
        ulong n;
        n = jacobi_one(G, chi1->conductor);
        if (chi1->parity)
            acb_set_si(res, -n);
        else
            acb_set_si(res, n);
    }
    else
    {
        if (G->q <= 150)
            acb_dirichlet_jacobi_sum_naive(res, G, chi1, chi2, prec);
        else if (G->num > 1)
            acb_dirichlet_jacobi_sum_primes(res, G, chi1, chi2, prec);
        else if (G->P[0].e > 1)
            acb_dirichlet_jacobi_sum_naive(res, G, chi1, chi2, prec);
        else
            acb_dirichlet_jacobi_sum_gauss(res, G, chi1, chi2, prec);
    }
}
