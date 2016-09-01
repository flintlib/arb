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

int main()
{
    slong prec = 64;
    n_primes_t iter;
    ulong q;

    flint_printf("thetanull....");
    fflush(stdout);

    /* look for vanishing theta values for prime power moduli */
    n_primes_init(iter);

#if 0 /* just to check q=300 & q=600 values */
    q = 2;
    while (q++ < 1000)
#else
    while ((q = n_primes_next(iter)) < 10000)
#endif
    {
        slong s;
        ulong k;
        acb_dirichlet_group_t G;
        acb_dirichlet_conrey_t x;
        acb_dirichlet_group_init(G, q);
        acb_dirichlet_conrey_init(x, G);
        acb_ptr theta, z;
        arb_t z1;
        arb_ptr t;

        theta = _acb_vec_init(G->phi_q);
        z =  _acb_vec_init(G->phi_q);

        arb_init(z1);
        arb_const_pi(z1, prec);
        arb_div_ui(z1, z1, q, prec);
        arb_neg(z1, z1);
        arb_exp(z1, z1, prec);

        t = _arb_vec_init(q);
        acb_dirichlet_arb_quadratic_powers(t, q, z1, prec);

        for (s = 0; s <= 1; s++)
        {
            k = 0;
            acb_dirichlet_conrey_one(x, G);
            do {
                acb_set_arb(z + k, t + x->n);
                k++;
            } while (acb_dirichlet_conrey_next(x, G) >= 0);

            acb_dirichlet_dft_conrey(theta, z, G, prec);

            for (k = 0; k < G->phi_q; k++)
            {
                if (acb_contains_zero(theta + k))
                {

                    slong i;
                    acb_dirichlet_conrey_one(x, G);
                    for (i = 0; i < k; i++)
                        acb_dirichlet_conrey_next(x, G);
                    if (acb_dirichlet_conrey_conductor(G,x) < q || acb_dirichlet_conrey_parity(G, x) != s)
                        continue;
                    flint_printf("\ntheta null q = %wu, n = %wu\n",q, x->n);
                    acb_printd(theta + k, 10);
                }
            }

            /* change t for odd characters */
            for (k = 0; k < q; k++)
                arb_mul_ui(t + k, t + k, k, prec);

        }

        _arb_vec_clear(t, q);
        _acb_vec_clear(theta, G->phi_q);
        _acb_vec_clear(z, G->phi_q);
        arb_clear(z1);
        acb_dirichlet_conrey_clear(x);
        acb_dirichlet_group_clear(G);
    }

    flint_cleanup();
    return EXIT_SUCCESS;
}
