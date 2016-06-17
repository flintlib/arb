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
        acb_dirichlet_group_t G;
        acb_dirichlet_conrey_t x;
        acb_dirichlet_char_t chi;
        ulong * v1, * v2, nv, k;

        acb_dirichlet_group_init(G, q);
        acb_dirichlet_conrey_init(x, G);
        acb_dirichlet_char_init(chi, G);

        nv = 100;
        v1 = flint_malloc(nv * sizeof(ulong));
        v2 = flint_malloc(nv * sizeof(ulong));

        acb_dirichlet_conrey_one(x, G);

        while (1) {

            acb_dirichlet_char_conrey(chi, G, x);

            acb_dirichlet_ui_chi_vec_loop(v1, G, chi, nv);
            acb_dirichlet_ui_chi_vec_primeloop(v2, G, chi, nv);

            if ((k = vec_diff(v1, v2, nv)))
            {
                flint_printf("FAIL: chi(%wu,%wu) = %wu != chi(%wu,%wu) = %wu mod %wu (modulus = %wu)\n",
                        chi->x->n, k, v1[k], chi->x->n, k, v2[k], chi->order, q);
                abort();
            }

            if (acb_dirichlet_conrey_next(x, G) == G->num)
                break;
        }

        flint_free(v1);
        flint_free(v2);
        acb_dirichlet_group_clear(G);
        acb_dirichlet_char_clear(chi);
        acb_dirichlet_conrey_clear(x);
    }

    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
