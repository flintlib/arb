/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

static ulong
vec_diff(ulong * v, ulong * ref, ulong nv)
{
    ulong k;
    for (k = 1; k < nv; k++)
        if (ref[k] != v[k])
            return k;
    return 0;
}

int main()
{
    ulong q;

    flint_printf("vec....");
    fflush(stdout);

    for (q = 2; q < 600; q ++)
    {
        dirichlet_group_t G;
        dirichlet_conrey_t x;
        dirichlet_char_t chi;
        ulong * v1, * v2, nv, k;

        dirichlet_group_init(G, q);
        dirichlet_conrey_init(x, G);
        dirichlet_char_init(chi, G);

        nv = 100;
        v1 = flint_malloc(nv * sizeof(ulong));
        v2 = flint_malloc(nv * sizeof(ulong));

        dirichlet_conrey_one(x, G);

        do {

            dirichlet_char_conrey(chi, G, x);

            dirichlet_ui_chi_vec_loop(v1, G, chi, nv);
            dirichlet_ui_chi_vec_primeloop(v2, G, chi, nv);

            if ((k = vec_diff(v1, v2, nv)))
            {
                flint_printf("FAIL: chi(%wu,%wu) = %wu != chi(%wu,%wu) = %wu mod %wu (modulus = %wu)\n",
                        chi->x->n, k, v1[k], chi->x->n, k, v2[k], chi->order, q);
                abort();
            }

        } while (dirichlet_conrey_next(x, G) >= 0);

        flint_free(v1);
        flint_free(v2);
        dirichlet_group_clear(G);
        dirichlet_char_clear(chi);
        dirichlet_conrey_clear(x);
    }

    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
