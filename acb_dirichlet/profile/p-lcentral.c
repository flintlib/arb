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

#define OUT 0

int main()
{
    slong prec = 100, digits = 30;
    ulong q;
    acb_t s;

    fflush(stdout);

    acb_init(s);
    acb_one(s);
    acb_div_si(s, s, 2, prec);

    /* look for vanishing theta values for prime power moduli */
    for (q = 3; q < 100; q++)
    {
        ulong k;
        acb_dirichlet_group_t G;
        acb_dirichlet_conrey_t x;
        acb_ptr z;

        if (q % 4 == 2)
            continue;

        acb_dirichlet_group_init(G, q);
        acb_dirichlet_conrey_init(x, G);

        z =  _acb_vec_init(G->phi_q);

        acb_dirichlet_l_vec_hurwitz(z, s, G, prec);

#if OUT
        k = 0;
        acb_dirichlet_conrey_one(x, G);
        while (acb_dirichlet_conrey_next(x, G) >= 0)
        { 
            k++;
            if (acb_dirichlet_conrey_conductor(G,x) < q)
                continue;
            flint_printf("%wu,%wu: ", q, x->n);
            acb_printd(z + k, digits);
            flint_printf("\n");
        }
#endif

        _acb_vec_clear(z, G->phi_q);
        acb_dirichlet_conrey_clear(x);
        acb_dirichlet_group_clear(G);
    }
    acb_clear(s);

    flint_cleanup();
    return EXIT_SUCCESS;
}
