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

    flint_printf("gauss....");
    fflush(stdout);

    /* check Gauss sums */

    for (q = 3; q < 250; q ++)
    {
        acb_dirichlet_group_t G;
        acb_dirichlet_conrey_t x;
        acb_dirichlet_char_t chi;

        acb_t s1, s2, s3, s4;

        acb_dirichlet_group_init(G, q);
        acb_dirichlet_conrey_init(x, G);
        acb_dirichlet_char_init(chi, G);

        acb_init(s1);
        acb_init(s2);
        acb_init(s3);
        acb_init(s4);
        acb_dirichlet_conrey_one(x, G);

        while (1) {

            acb_dirichlet_char_conrey(chi, G, x);

            acb_dirichlet_gauss_sum_naive(s1, G, chi, prec);
            acb_dirichlet_gauss_sum(s2, G, chi, prec);
            acb_dirichlet_gauss_sum_factor(s3, G, chi, prec);
            if (chi->conductor == G->q)
                acb_dirichlet_gauss_sum_theta(s4, G, chi, prec);
            else
                acb_set(s4, s1);

            if (!acb_overlaps(s1, s2)
                    || !acb_overlaps(s1, s3)
                    || !acb_overlaps(s1, s4))
            {
                flint_printf("FAIL: G(chi_%wu(%wu))\n\n", q, chi->x->n);
                flint_printf("\nnaive ", q, x->n);
                acb_printd(s1, 25);
                flint_printf("\ndefault ", q, x->n);
                acb_printd(s2, 25);
                flint_printf("\nfactor ", q, x->n);
                acb_printd(s3, 25);
                flint_printf("\ntheta ", q, x->n);
                acb_printd(s4, 25);
                abort();
            }

            if (acb_dirichlet_conrey_next(x, G) == G->num)
                break;
        }
        acb_clear(s1);
        acb_clear(s2);
        acb_clear(s3);
        acb_clear(s4);

        acb_dirichlet_group_clear(G);
        acb_dirichlet_char_clear(chi);
        acb_dirichlet_conrey_clear(x);
    }

    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
