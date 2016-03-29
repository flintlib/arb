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
#include <math.h>

#define PI   3.14159265358
#define LOG2 0.69314718055

int main()
{
    slong prec = 80;
    ulong q;

    flint_printf("thetanull....");
    fflush(stdout);

    /* check the only theta functions
     *   theta(chi) = sum chi(k)* k^odd * exp(-Pi * k^2 / q)
     * vanishing at 1 correspond to two specific
     * characters of moduli 300 and 600 + conjugates
     */

    for (q = 3; q < 800; q ++)
    {
        acb_dirichlet_group_t G;
        acb_dirichlet_conrey_t x;
        acb_dirichlet_char_t chi;
        ulong * v, nv, k;

        acb_t zeta, sum;
        acb_ptr z;

        arb_t eq;
        arb_ptr t, kt, tt;

        if (q % 4 == 2)
            continue;

        acb_dirichlet_group_init(G, q);
        acb_dirichlet_conrey_init(x, G);
        acb_dirichlet_char_init(chi, G);

        acb_init(zeta);
        acb_dirichlet_zeta(zeta, G->expo, prec);

        z = _acb_vec_init(G->expo);
        _acb_vec_set_powers(z, zeta, G->expo, prec);

        nv = ceil(sqrt((double)q * prec * LOG2 / PI)) + 2;
        v = flint_malloc(nv * sizeof(ulong));

        arb_init(eq);
        arb_const_pi(eq, prec);
        arb_div_ui(eq, eq, q, prec);
        arb_neg(eq, eq);
        arb_exp(eq, eq, prec);

        t = _arb_vec_init(nv);
        acb_dirichlet_arb_quadratic_powers(t, nv, eq, prec);

        kt = _arb_vec_init(nv);
        for (k = 1; k < nv; k++)
            arb_mul_ui(kt + k, t + k, k, prec);

        /* theta function on primitive characters */
        acb_init(sum);
        acb_dirichlet_conrey_first_primitive(x, G);

        while (1) {
            ulong m;
            acb_zero(sum);

            acb_dirichlet_char_conrey(chi, G, x);
            acb_dirichlet_chi_vec_loop(v, nv, G, chi);

            m = G->expo / chi->order;
            /*
            flint_printf("Theta(chi_%wu(%wu)) (m=%wu)\n", q, chi->n, m);
            */
            tt = acb_dirichlet_char_parity(chi) ? kt : t;

            for (k = 1; k < nv; k++)
                if (v[k] != ACB_DIRICHLET_CHI_NULL)
                    acb_addmul_arb(sum, z + (v[k] * m), tt + k, prec);

            if ((q == 300 && (chi->n == 271 || chi->n == 131))
                    || (q == 600 && (chi->n == 11 || chi->n == 91)))
            {
                if (!acb_contains_zero(sum))
                {
                flint_printf("FAIL: Theta(chi_%wu(%wu))=", q, chi->n);
                acb_printd(sum, 10);
                flint_printf("\n");
                abort();
                }
            }
            else if (acb_contains_zero(sum))
            {
                flint_printf("FAIL: Theta(chi_%wu(%wu))=", q, chi->n);
                acb_printd(sum, 10);
                flint_printf("\n");
                abort();
            }

            if (acb_dirichlet_conrey_next_primitive(x, G) == G->num)
                break;
        }
        _acb_vec_clear(z, G->expo);
        _arb_vec_clear(t, nv);
        acb_clear(zeta);
        acb_clear(sum);
        arb_clear(eq);
        acb_dirichlet_group_clear(G);
        acb_dirichlet_char_clear(chi);
        acb_dirichlet_conrey_clear(x);
    }

    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
