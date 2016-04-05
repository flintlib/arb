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
    slong iter;
    flint_rand_t state;

    flint_printf("conrey....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 3000; iter++)
    {
        acb_dirichlet_group_t G;
        acb_dirichlet_conrey_t x;
        ulong q, n, k, sum;
        long ref;
        /*int * bits;*/

        q = 1 + n_randint(state, 1000 * (1 + iter / 100));

        acb_dirichlet_group_init(G, q);

        acb_dirichlet_conrey_init(x, G);

        /* check group size and elements */
        acb_dirichlet_conrey_one(x, G);
        sum = 1;

#if 1
        for (n = 1; acb_dirichlet_conrey_next(x, G) < G->num; n++)
            sum += x->n * x->n;
#else
        /* iteration much faster than gcd below */
        n = 1;
        for (k = 2; k < G->q; k++)
        {
            if (n_gcd(k, G->q) > 1)
                continue;
            n++;
            sum += k * k;
        }
#endif

        /* use http://oeis.org/A053818 to check all elements
         * are gone through */
        ref = (q % 4 == 2) ? -2 : 1;
        k = (G->neven == 2) ? 1 : 0;
        for (; k<G->num; k++)
          ref = - ref * G->primes[k];
        ref = ( G->phi_q * (2 * q * q + ref) ) / 6;


        if (n != G->phi_q)
        {
            flint_printf("FAIL: group size\n\n");
            flint_printf("q = %wu\n\n", q);
            flint_printf("phi(q) = %wu\n\n", G->phi_q);
            flint_printf("loop index = %wu\n\n", n);
            abort();
        }
        if (sum != ref && q > 1)
        {
            flint_printf("FAIL: sum test\n\n");
            flint_printf("q = %wu\n\n", q);
            flint_printf("sum k^2 = %wu\n\n", ref);
            flint_printf("sum obtained = %wu\n\n", sum);
            abort();
        }

        if (q % 4 == 2)
            continue;

        /* check primitive elements */
        if (q % 4 == 2)
            ref = 0;
        else
        {
           ref = 1;
           k = (G->neven == 2) ? 1 : 0;
           for (; k < G->num; k++)
           {
               ulong p = G->primes[k];
               if (G->exponents[k] == 1)
                   ref *= p - 2;
               else
                   ref *= (p * (p - 2) + 1) * n_pow(p, G->exponents[k] - 2);
           }
        }

        acb_dirichlet_conrey_first_primitive(x, G);
        for (n = 1; (k=acb_dirichlet_conrey_next_primitive(x, G)) < G->num; n++);

        if (n != ref)
        {
            flint_printf("FAIL: number of primitive elements\n\n");
            flint_printf("q = %wu\n\n", q);
            flint_printf("# primitive = %wu\n\n", ref);
            flint_printf("loop index = %wu\n\n", n);
            abort();
        }

        acb_dirichlet_conrey_clear(x);
        acb_dirichlet_group_clear(G);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
