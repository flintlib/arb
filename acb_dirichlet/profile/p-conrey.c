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
#include "profiler.h"

int main()
{
    slong iter, k, n, nref;
    ulong q, maxq;

    maxq = 4000;
    flint_printf("loop over all (Z/q)* for q<=%wu....\n", maxq);
    fflush(stdout);

    flint_printf("gcd................ ");
    TIMEIT_ONCE_START
    for (n = 0, q = 2; q <= maxq; q++)
        for (k = 1; k < q; k++)
            if (n_gcd(k, q) > 1)
                continue;
            else
                n++;
    TIMEIT_ONCE_STOP
    nref = n;

    flint_printf("conrey.... ");
    TIMEIT_ONCE_START
    for (n = 0, q = 2; q <= maxq; q++)
    {
        acb_dirichlet_group_t G;
        acb_dirichlet_conrey_t x;

        acb_dirichlet_group_init(G, q);
        acb_dirichlet_conrey_init(x, G);

        acb_dirichlet_conrey_one(x, G);
        n++;

        for (; acb_dirichlet_conrey_next(x, G) >= 0; n++);
        acb_dirichlet_conrey_clear(x);
        acb_dirichlet_group_clear(G);
    }
    TIMEIT_ONCE_STOP
    if (n != nref)
    {
            flint_printf("FAIL: wrong number of elements %wu != %wu\n\n",n, nref);
            abort();
    }

    flint_printf("chars.... ");
    TIMEIT_ONCE_START
    for (n = 0, q = 2; q <= maxq; q++)
    {
        acb_dirichlet_group_t G;
        acb_dirichlet_char_t chi;

        acb_dirichlet_group_init(G, q);
        acb_dirichlet_char_init(chi, G);

        acb_dirichlet_char_one(chi, G);
        n++;

        for (; acb_dirichlet_char_next(chi, G) >= 0; n++);
        acb_dirichlet_char_clear(chi);
        acb_dirichlet_group_clear(G);
    }
    TIMEIT_ONCE_STOP
    if (n != nref)
    {
            flint_printf("FAIL: wrong number of elements %wu != %wu\n\n",n, nref);
            abort();
    }


    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
