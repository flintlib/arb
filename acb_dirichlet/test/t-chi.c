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
        acb_t zn1, zn2, zn1n2, zn1zn2;
        acb_dirichlet_group_t G;
        ulong q, m, n1, n2, iter2;

        q = 1 + n_randint(state, 1000);

        acb_dirichlet_group_init(G, q);
        acb_init(zn1);
        acb_init(zn2);
        acb_init(zn1n2);
        acb_init(zn1zn2);

        for (iter2 = 0; iter2 < 10; iter2++)
        {
            do {
                m = 1 + n_randint(state, q);
            } while (n_gcd(q, m) != 1);

            n1 = n_randint(state, 1000);
            n2 = n_randint(state, 1000);

            acb_dirichlet_chi(zn1, G, m, n1, 53);
            acb_dirichlet_chi(zn2, G, m, n2, 53);
            acb_dirichlet_chi(zn1n2, G, m, n1 * n2, 53);
            acb_mul(zn1zn2, zn1, zn2, 53);

            if (!acb_overlaps(zn1n2, zn1zn2))
            {
                flint_printf("FAIL: overlap\n\n");
                flint_printf("q = %wu\n\n", q);
                flint_printf("m = %wu\n\n", m);
                flint_printf("n1 = %wu\n\n", n1);
                flint_printf("n2 = %wu\n\n", n2);
                flint_printf("zn1 = "); acb_printd(zn1, 15); flint_printf("\n\n");
                flint_printf("zn2 = "); acb_printd(zn2, 15); flint_printf("\n\n");
                flint_printf("zn1n2 = "); acb_printd(zn1n2, 15); flint_printf("\n\n");
                flint_printf("zn1zn2 = "); acb_printd(zn1zn2, 15); flint_printf("\n\n");
                abort();
            }
        }

        acb_dirichlet_group_clear(G);
        acb_clear(zn1);
        acb_clear(zn2);
        acb_clear(zn1n2);
        acb_clear(zn1zn2);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

