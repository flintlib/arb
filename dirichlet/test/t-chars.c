/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "dirichlet.h"

int main()
{
    slong iter, bits;
    flint_rand_t state;

    flint_printf("chars....");
    fflush(stdout);
    flint_randinit(state);
    for (bits = 5; bits <= 30; bits += 5)
    {

        for (iter = 0; iter < 50; iter++)
        {
            dirichlet_group_t G;
            dirichlet_conrey_t x;
            dirichlet_char_t chi, chi2;
            ulong q, iter2;

            q = 2 + n_randint(state, 1 << bits);

            dirichlet_group_init(G, q);
            dirichlet_conrey_init(x, G);
            dirichlet_char_init(chi, G);
            dirichlet_char_init(chi2, G);

            dirichlet_group_dlog_precompute(G, 50);

            /* check number char properties */
            for (iter2 = 0; iter2 < 100; iter2++)
            {
                int par;
                ulong m, n;
                ulong order, chim1, pairing, cn, cm, cond;

                do
                    m = n_randint(state, q);
                while (n_gcd(q, m) > 1);

                dirichlet_char(chi, G, m);
                dirichlet_conrey_log(x, G, m);
                dirichlet_char_conrey(chi2, G, x);

                if (!dirichlet_char_eq_deep(G, chi, chi2))
                {
                    flint_printf("FAIL: init char\n\n");
                    flint_printf("q = %wu\n\n", q);
                    flint_printf("m = %wu\n\n", m);
                    dirichlet_char_print(G, chi);
                    flint_printf("\n");
                    dirichlet_char_print(G, chi2);
                    flint_printf("\n");
                    abort();
                }

                order = dirichlet_ui_order(G, m);
                if (order != chi->order.n)
                {
                    flint_printf("FAIL: order\n\n");
                    flint_printf("q = %wu\n\n", q);
                    flint_printf("m = %wu\n\n", m);
                    flint_printf("order(m) = %wu\n\n", order);
                    flint_printf("chi->order = %wu\n\n", chi->order);
                    abort();
                }

                cond = dirichlet_ui_conductor(G, m);
                if (cond != chi->conductor)
                {
                    flint_printf("FAIL: conductor\n\n");
                    flint_printf("q = %wu\n\n", q);
                    flint_printf("m = %wu\n\n", m);
                    flint_printf("conductor(m) = %wu\n\n", cond);
                    flint_printf("chi->conductor = %wu\n\n", chi->conductor);
                    abort();
                }

                par = dirichlet_ui_parity(G, m);
                chim1 = dirichlet_ui_chi(G, chi, q - 1);
                if (dirichlet_char_parity(chi) != par || par != (chim1 != 0))
                {
                    flint_printf("FAIL: parity\n\n");
                    flint_printf("q = %wu\n\n", q);
                    flint_printf("m = %wu\n\n", m);
                    flint_printf("chi(-1) = %wu\n\n", chim1);
                    flint_printf("char_parity = %d", dirichlet_char_parity(chi));
                    flint_printf("parity_ui = %d", par);
                    dirichlet_char_print(G, chi);
                    abort();
                }

                do
                    n = n_randint(state, q);
                while (n_gcd(q, n) > 1);

                dirichlet_char(chi2, G, n);
                pairing = dirichlet_ui_pairing(G, m, n);
                cn = dirichlet_ui_chi(G, chi, n) * (G->expo / chi->order.n);
                cm = dirichlet_ui_chi(G, chi2, m) * (G->expo / chi2->order.n);

                if (pairing != cn || pairing != cm)
                {
                    flint_printf("FAIL: pairing\n\n");
                    flint_printf("q = %wu\n\n", q);
                    flint_printf("m = %wu\n\n", m);
                    flint_printf("n = %wu\n\n", n);
                    flint_printf("chi(m,n) = %wu\n\n", pairing);
                    flint_printf("chi(m)(n) = %wu\n\n", cn);
                    flint_printf("chi(n)(m) = %wu\n\n", cm);
                    abort();
                }

                dirichlet_conrey_next(x, G);
                dirichlet_char_next(chi, G);
                dirichlet_char_conrey(chi2, G, x);

                if (!dirichlet_char_eq_deep(G, chi, chi2))
                {
                    flint_printf("FAIL: next char\n\n");
                    flint_printf("q = %wu\n\n", q);
                    flint_printf("m = %wu\n\n", m);
                    dirichlet_char_print(G, chi);
                    flint_printf("\n");
                    dirichlet_char_print(G, chi2);
                    flint_printf("\n");
                    abort();
                }

            }

            dirichlet_group_dlog_clear(G);

            dirichlet_char_clear(chi);
            dirichlet_char_clear(chi2);
            dirichlet_conrey_clear(x);
            dirichlet_group_clear(G);
        }
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
