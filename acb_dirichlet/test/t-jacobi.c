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
    slong prec = 128;
    ulong q;

    flint_printf("jacobi....");
    fflush(stdout);

    /* check Jacobi sums */

    for (q = 3; q < 250; q ++)
    {
        acb_dirichlet_group_t G;
        acb_dirichlet_char_t chi1, chi2;

        acb_t s1, s2;

        acb_dirichlet_group_init(G, q);
        acb_dirichlet_char_init(chi1, G);
        acb_dirichlet_char_init(chi2, G);

        acb_init(s1);
        acb_init(s2);

        acb_dirichlet_char_one(chi1, G);

        while (1) {

            acb_dirichlet_char_one(chi2, G);

            while (1) {

                acb_dirichlet_jacobi_sum_naive(s1, G, chi1, chi2, prec);
                acb_dirichlet_jacobi_sum(s2, G, chi1, chi2, prec);

                if (!acb_overlaps(s1, s2))
                {
                    flint_printf("FAIL: J_%wu(%wu,%wu)",
                            q, chi1->x->n, chi2->x->n);
                    flint_printf("\nnaive ");
                    acb_printd(s1, 25);
                    flint_printf("\ntheta ");
                    acb_printd(s2, 25);
                    flint_printf("\n");
                    flint_printf("cond = %wu, %wu, %wu\n",
                            chi1->conductor, chi2->conductor,
                            acb_dirichlet_ui_conductor(G, nmod_mul(chi1->x->n, chi2->x->n, G->mod))
                            );
                    abort();
                }
                if (acb_dirichlet_char_next(chi2, G) == G->num)
                    break;

            }

            if (acb_dirichlet_char_next(chi1, G) == G->num)
                break;

        }

        acb_clear(s1);
        acb_clear(s2);
        acb_dirichlet_group_clear(G);
        acb_dirichlet_char_clear(chi1);
        acb_dirichlet_char_clear(chi2);
    }

    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
