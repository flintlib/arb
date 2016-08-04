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
        acb_dirichlet_conrey_t x, y;
        ulong q, n, k, sum;
        long ref;
        /*int * bits;*/

        q = 1 + n_randint(state, 1000 * (1 + iter / 100));

        acb_dirichlet_group_init(G, q);

        acb_dirichlet_conrey_init(x, G);
        acb_dirichlet_conrey_init(y, G);

        /* check group size and elements */
        acb_dirichlet_conrey_one(x, G);
        sum = 1;

#if 1
        for (n = 1; acb_dirichlet_conrey_next(x, G) >= 0; n++)
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
        for (k = (G->neven == 2); k < G->num; k++)
            ref = - ref * G->P[k].p;
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

        acb_dirichlet_conrey_first_primitive(x, G);
        for (n = 1; acb_dirichlet_conrey_next_primitive(x, G) >= 0; n++);

        ref = acb_dirichlet_number_primitive(G);
        if (n != ref)
        {
            flint_printf("FAIL: number of primitive elements\n\n");
            flint_printf("q = %wu\n\n", q);
            flint_printf("# primitive = %wu\n\n", ref);
            flint_printf("loop index = %wu\n\n", n);
            abort();
        }


        /* some random elements, check log and exp */
        for (n = 0; n < 30; n++)
        {
            slong k;
            ulong m;
            for (m = 1; n_gcd(m, q) > 1; m = n_randint(state, q));
            acb_dirichlet_conrey_log(x, G, m);

            if (m != acb_dirichlet_conrey_exp(x, G))
            {
                flint_printf("FAIL: conrey log and exp\n\n");
                flint_printf("q = %wu\n\n", q);
                flint_printf("m = %wu\n\n", m);
                flint_printf("conrey = ");
                acb_dirichlet_conrey_print(G, x);
                flint_printf("\n\nnumber = %wu\n\n", x->n);
                abort();
            }

            for (k = 0; k < G->num; k++)
                x->log[k] = n_randint(state, G->P[k].phi);

            m = acb_dirichlet_conrey_exp(x, G);
            acb_dirichlet_conrey_log(y, G, m);

            if (!acb_dirichlet_conrey_eq(G, x, y))
            {
                flint_printf("FAIL: conrey log and exp\n\n");
                flint_printf("q = %wu\n\n", q);
                flint_printf("conrey = ");
                acb_dirichlet_conrey_print(G, x);
                flint_printf("m = %wu\n\n", m);
                flint_printf("log = ");
                acb_dirichlet_conrey_print(G, y);
                flint_printf("\n\nnumber = %wu\n\n", y->n);
                abort();
            }

        }

        acb_dirichlet_conrey_clear(x);
        acb_dirichlet_conrey_clear(y);
        acb_dirichlet_group_clear(G);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
