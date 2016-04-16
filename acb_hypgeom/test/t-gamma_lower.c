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

    Copyright (C) 2016 Arb authors

******************************************************************************/

#include "acb_hypgeom.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("gamma_lower....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 2000 * arb_test_multiplier(); iter++)
    {
        acb_t a0, a1, b, z, w0, w1, t, u;
        slong prec0, prec1;
        int regularized;

        acb_init(a0);
        acb_init(a1);
        acb_init(b);
        acb_init(z);
        acb_init(w0);
        acb_init(w1);
        acb_init(t);
        acb_init(u);

        regularized = n_randint(state, 3);

        prec0 = 2 + n_randint(state, 1000);
        prec1 = 2 + n_randint(state, 1000);

        acb_randtest_param(a0, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(z, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(w0, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(w1, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        acb_add_ui(a1, a0, 1, prec0);

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_gamma_lower_1f1a(w0, a0, z, regularized, prec0);
                break;
            case 1:
                acb_hypgeom_gamma_lower_1f1b(w0, a0, z, regularized, prec0);
                break;
            default:
                acb_hypgeom_gamma_lower(w0, a0, z, regularized, prec0);
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_gamma_lower_1f1a(w1, a0, z, regularized, prec1);
                break;
            case 1:
                acb_hypgeom_gamma_lower_1f1b(w1, a0, z, regularized, prec1);
                break;
            default:
                acb_hypgeom_gamma_lower(w1, a0, z, regularized, prec1);
        }

        if (!acb_overlaps(w0, w1))
        {
            flint_printf("FAIL: consistency\n\n");
            flint_printf("a0 = "); acb_printd(a0, 30); flint_printf("\n\n");
            flint_printf("z = "); acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("w0 = "); acb_printd(w0, 30); flint_printf("\n\n");
            flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
            abort();
        }

        switch (n_randint(state, 3))
        {
            case 1:
                acb_hypgeom_gamma_lower_1f1a(w1, a1, z, regularized, prec1);
                break;
            case 2:
                acb_hypgeom_gamma_lower_1f1b(w1, a1, z, regularized, prec1);
                break;
            default:
                acb_hypgeom_gamma_lower(w1, a1, z, regularized, prec1);
        }

        if (regularized == 2)
        {
            /* a r(a,z) - exp(-z) - z r(a+1,z) = 0 */
            acb_neg(u, z);
            acb_exp(u, u, prec0);
            acb_neg(t, u);

            acb_mul(b, w1, z, prec0);
            acb_addmul(t, a0, w0, prec0);
            acb_sub(t, t, b, prec0);
        }
        else if (regularized == 1)
        {
            /* q(a,z) - exp(-z) z^a / Gamma(a+1) - q(a+1,z) = 0 */
            acb_pow(t, z, a0, prec0);
            acb_rgamma(u, a1, prec0);
            acb_mul(t, t, u, prec0);

            acb_neg(u, z);
            acb_exp(u, u, prec0);
            acb_neg(u, u);
            acb_mul(t, t, u, prec0);

            acb_add(t, t, w0, prec0);
            acb_sub(t, t, w1, prec0);
        }
        else
        {
            /* a gamma(a,z) - exp(-z) z^a - gamma(a+1,z) = 0 */
            acb_pow(t, z, a0, prec0);

            acb_neg(u, z);
            acb_exp(u, u, prec0);
            acb_neg(u, u);
            acb_mul(t, t, u, prec0);

            acb_addmul(t, a0, w0, prec0);
            acb_sub(t, t, w1, prec0);
        }

        if (!acb_contains_zero(t))
        {
            flint_printf("FAIL: contiguous relation\n\n");
            flint_printf("regularized = %d\n\n", regularized);
            flint_printf("a0 = "); acb_printd(a0, 30); flint_printf("\n\n");
            flint_printf("z = ");  acb_printd(z, 30); flint_printf("\n\n");
            flint_printf("w0 = "); acb_printd(w0, 30); flint_printf("\n\n");
            flint_printf("w1 = "); acb_printd(w1, 30); flint_printf("\n\n");
            flint_printf("t = "); acb_printd(t, 30); flint_printf("\n\n");
            abort();
        }

        acb_clear(a0);
        acb_clear(a1);
        acb_clear(b);
        acb_clear(z);
        acb_clear(w0);
        acb_clear(w1);
        acb_clear(t);
        acb_clear(u);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
