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

    Copyright (C) 2016 Fredrik Johansson

******************************************************************************/

#include "acb_dirichlet.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("chi....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        acb_dirichlet_group_t G;
        acb_conrey_t x;
        ulong q, n, k, sum;
        long ref;
        /*int * bits;*/

        q = 1 + n_randint(state, 1000 * (1 + iter / 100));

        acb_dirichlet_group_init(G, q);

        /* use http://oeis.org/A053818 to check all elements
         * are gone through */
        ref = (q % 4 == 2) ? -2 : 1;
        k = (G->neven == 2) ? 1 : 0;
        for (; k<G->num; k++)
          ref = - ref * G->primes[k];
        ref = ( G->phi_q * (2 * q * q + ref) ) / 6;

        acb_conrey_init(x, G);
        acb_conrey_one(x, G);
        sum = 1;

        /* check group size */
        for (n = 1; acb_conrey_next(x, G) < G->num; n++)
            sum += x->n * x->n;

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

        acb_conrey_clear(x);
        acb_dirichlet_group_clear(G);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
