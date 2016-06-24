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
charsum_1modsomep(const acb_dirichlet_group_t G, ulong cond)
{
    slong k, f = 1, mu = 1, pow = 1;

    for (k = 0; k < G->num; k++)
    {
        ulong p = G->primes[k];
        if (G->exponents[k] > 1)
        {
            if (cond % (p*p))
                pow *= G->primepowers[k] / p;
            else
                return 0;
        }
        if (cond % p == 0)  /* p | conductor */
            mu *= -1;
        else
            f *= p - 2;
    }

    return mu * pow * f;
}

void
acb_dirichlet_jacobi_sum(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi1, const acb_dirichlet_char_t chi2, slong prec)
{

    if (G->q_even > 1)
    {
        acb_zero(res);
        return;
    }

    if (chi1->x->n == 1 || chi2->x->n == 1)
    {
        if (chi1->x->n == 1 && chi2->x->n == 1)
        {
            /* q = prod p^e -> prod (p-2)p^(e-1) invertible x & 1-x */
            slong k, n = 1;
            //flint_printf("## a=b=1[q]\n");
            for (k = 0; k < G->num; k++)
            {
                n *= G->primes[k] - 2;
                if (G->exponents[k] > 1)
                    n *= G->primepowers[k] / G->primes[k];
            }
            acb_set_si(res, n);
        }
        else
        {
            ulong cond = (chi1->x->n == 1) ? chi2->conductor : chi1->conductor;
            acb_set_si(res, charsum_1modsomep(G, cond));
        }
    }
    else if (nmod_mul(chi1->x->n, chi2->x->n, G->mod) == 1)
    {
        ulong n;
        //flint_printf("## ab=1[q]\n");
        n = charsum_1modsomep(G, chi1->conductor);
        if (chi1->parity)
            acb_set_si(res, -n);
        else
            acb_set_si(res, n);
    }
    else
    {
        /* J_q(a,b)G_q(ab) = G_q(a)G_q(b) */
        acb_dirichlet_char_t chi12;

        //flint_printf("## via gauss\n");

        acb_dirichlet_char_init(chi12, G);
        acb_dirichlet_char_mul(chi12, G, chi1, chi2);

        if (chi12->conductor != G->q)
        {
            //flint_printf("## jacobi: non primitive product, %wu * %wu -> %wu\n",
                    //chi1->conductor, chi2->conductor, chi12->conductor);
            //acb_dirichlet_jacobi_sum_naive(res, G, chi1, chi2, prec);
        }

        if (1)
        {
            acb_t tmp;
            acb_init(tmp);

            /* FIXME: remove naive */
            acb_dirichlet_gauss_sum_naive(res, G, chi1, prec);
            acb_dirichlet_gauss_sum_naive(tmp, G, chi2, prec);
            acb_mul(res, res, tmp, prec);
            acb_dirichlet_gauss_sum_naive(tmp, G, chi12, prec);
            acb_div(res, res, tmp, prec);
            if (chi12->conductor < G->q)
            {
                /* à la louche... */
                if (chi1->conductor == chi2->conductor
                        && chi2->conductor == chi12->conductor)
                {
                    slong k;
                    slong m = 1;
                    for (k = 0; k < G->num; k++)
                    {
                        ulong p = G->primes[k];
                        if (chi1->conductor % p) 
                            m = - m * (p - 2);
                    }
                    /*
                    flint_printf("cond = %wu, %wu, %wu -> mult by %wd\n",
                            chi1->conductor, chi2->conductor, chi12->conductor,
                            m);
                    */
                    acb_mul_si(res, res, m, prec); 
                }
                else
                    acb_div_si(res, res, G->q / chi12->conductor, prec);
            }


            acb_dirichlet_char_clear(chi12);
            acb_clear(tmp);
        }
    }

}
