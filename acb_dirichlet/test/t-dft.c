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

    slong k;
    slong prec = 100;
    slong nq = 10;
    ulong q[10] = { 2, 3, 4, 5, 6, 10, 15, 30, 308, 961};

    flint_printf("dft....");
    fflush(stdout);

    /* cyclic dft */
    for (k = 0; k < nq; k++)
    {
        slong i;
        acb_ptr v, w1, w2;

        v = _acb_vec_init(q[k]);
        w1 = _acb_vec_init(q[k]);
        w2 = _acb_vec_init(q[k]);

        for (i = 0; i < q[k]; i++)
            acb_set_si(v + i, i);

        acb_dirichlet_dft_pol(w1, v, q[k], prec);
        acb_dirichlet_dft_fast(w2, v, q[k], prec);

        for (i = 0; i < q[k]; i++)
        {
            if (!acb_overlaps(w1 + i, w2 + i))
            {
                flint_printf("differ from index %ld / %ld \n\n",i,q[k]);
                flint_printf("pol =\n");
                acb_vec_printd(w1, q[k], 10);
                flint_printf("fast =\n");
                acb_vec_printd(w2, q[k], 10);
                flint_printf("\n\n");
                abort();
            }
        }

        _acb_vec_clear(v, q[k]);
        _acb_vec_clear(w1, q[k]);
        _acb_vec_clear(w2, q[k]);
    }

    /* Dirichlet group DFT */
    for (k = 0; k < nq - 1; k++)
    {
        slong i, j, len;
        acb_dirichlet_group_t G;
        acb_dirichlet_conrey_t x, y;
        acb_t chiy;
        acb_ptr v, w1, w2;

        acb_dirichlet_group_init(G, q[k]);

        len = G->phi_q;
        v = _acb_vec_init(len);
        w1 = _acb_vec_init(len);
        w2 = _acb_vec_init(len);

        acb_dirichlet_conrey_init(x, G);
        acb_dirichlet_conrey_one(x, G);
        for (i = 0; i < len; i++)
        {
            acb_set_si(v + i, x->n);
            acb_dirichlet_conrey_next(x, G);
        }

        /* naive */
        acb_init(chiy);
        acb_dirichlet_conrey_init(y, G);
        acb_dirichlet_conrey_one(x, G);
        for (i = 0;  i < len; i++)
        {
            acb_zero(w1 + i);
            acb_dirichlet_conrey_one(y, G);
            for (j = 0; j < len; j++)
            {
                acb_dirichlet_pairing_conrey(chiy, G, x, y, prec);
                acb_addmul(w1 + i, chiy, v + j, prec);
                acb_dirichlet_conrey_next(y, G);
            }
            acb_dirichlet_conrey_next(x, G);
        }
        acb_clear(chiy);
        acb_dirichlet_conrey_clear(y);
        acb_dirichlet_conrey_clear(x);

        /* dft */
        acb_dirichlet_dft_conrey(w2, v, G, prec);

        for (i = 0; i < len; i++)
        {
            if (!acb_overlaps(w1 + i, w2 + i))
            {
                flint_printf("FAIL\n\n");
                flint_printf("q = %wu\n", q[k]);
                flint_printf("v [size %wu]\n", len);
                acb_vec_printd(v, len, 10);
                flint_printf("\nDFT differ from index %ld / %ld \n", i, len);
                flint_printf("\nnaive =\n");
                acb_vec_printd(w1, len, 10);
                flint_printf("\nfast =\n");
                acb_vec_printd(w2, len, 10);
                flint_printf("\n\n");
                abort();
            }
        }

        _acb_vec_clear(v, len);
        _acb_vec_clear(w1, len);
        _acb_vec_clear(w2, len);

        acb_dirichlet_group_clear(G);
    }

    flint_printf("PASS\n");
}
