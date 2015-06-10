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

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include "acb_hypgeom.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("bessel_y....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 2000; iter++)
    {
        acb_t nu, nu1, z, jv, jv1, yv, yv1, r, s;
        long prec, n;

        acb_init(nu);
        acb_init(nu1);
        acb_init(z);
        acb_init(jv);
        acb_init(jv1);
        acb_init(yv);
        acb_init(yv1);
        acb_init(r);
        acb_init(s);

        prec = 2 + n_randint(state, 500);

        acb_randtest_param(nu, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 10));
        acb_randtest(z, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 10));

        if (n_randint(state, 4) == 0)
            arb_zero(acb_imagref(z));

        acb_add_ui(nu1, nu, 1, prec);

        acb_hypgeom_bessel_j(jv, nu, z, prec);
        acb_hypgeom_bessel_j(jv1, nu1, z, prec);
        acb_hypgeom_bessel_y(yv, nu, z, prec);
        acb_hypgeom_bessel_y(yv1, nu1, z, prec);

        acb_mul(r, jv1, yv, prec);
        acb_submul(r, jv, yv1, prec);
        acb_mul(r, r, z, prec);
        acb_const_pi(s, prec);
        acb_mul(r, r, s, prec);
        acb_sub_ui(r, r, 2, prec);

        if (!acb_contains_zero(r))
        {
            printf("FAIL: wronskian\n\n");
            printf("nu = "); acb_printd(nu, 30); printf("\n\n");
            printf("z = "); acb_printd(z, 30); printf("\n\n");
            printf("yv = "); acb_printd(yv, 30); printf("\n\n");
            printf("yv1 = "); acb_printd(yv1, 30); printf("\n\n");
            printf("jv = "); acb_printd(jv, 30); printf("\n\n");
            printf("jv1 = "); acb_printd(jv1, 30); printf("\n\n");
            printf("r = "); acb_printd(r, 30); printf("\n\n");
            abort();
        }

        /* Y_n(-z) = (-1)^n [Y_n(z) - (2/pi) [log(z) - log(-z)] J_v(z)] */
        n = n_randint(state, 20) - 10;
        acb_set_si(nu, n);

        acb_hypgeom_bessel_y(yv, nu, z, prec);
        acb_hypgeom_bessel_j(jv, nu, z, prec);

        acb_log(r, z, prec);
        acb_neg(s, z);
        acb_log(s, s, prec);
        acb_sub(r, r, s, prec);
        acb_mul(jv, jv, r, prec);
        acb_const_pi(r, prec);
        acb_div(jv, jv, r, prec);
        acb_mul_2exp_si(jv, jv, 1);

        acb_sub(r, yv, jv, prec);
        if (n % 2)
            acb_neg(r, r);

        acb_neg(yv1, z);
        acb_hypgeom_bessel_y(yv1, nu, yv1, prec);

        if (!acb_overlaps(r, yv1))
        {
            printf("FAIL: reflection formula\n\n");
            printf("nu = "); acb_printd(nu, 30); printf("\n\n");
            printf("z = "); acb_printd(z, 30); printf("\n\n");
            printf("yv = "); acb_printd(yv, 30); printf("\n\n");
            printf("yv1 = "); acb_printd(yv1, 30); printf("\n\n");
            printf("jv = "); acb_printd(jv, 30); printf("\n\n");
            printf("r = "); acb_printd(r, 30); printf("\n\n");
            abort();
        }

        acb_clear(nu);
        acb_clear(nu1);
        acb_clear(z);
        acb_clear(jv);
        acb_clear(jv1);
        acb_clear(yv);
        acb_clear(yv1);
        acb_clear(r);
        acb_clear(s);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

