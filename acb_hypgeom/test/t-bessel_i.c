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
    slong iter;
    flint_rand_t state;

    printf("bessel_i....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 2000; iter++)
    {
        acb_t nu, z, jv, iv, t;
        slong prec;

        acb_init(nu);
        acb_init(z);
        acb_init(jv);
        acb_init(iv);
        acb_init(t);

        prec = 2 + n_randint(state, 500);

        acb_randtest_param(nu, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 10));
        acb_randtest(z, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_bessel_i_asymp(iv, nu, z, prec);
                break;
            case 1:
                acb_hypgeom_bessel_i_0f1(iv, nu, z, prec);
                break;
            default:
                acb_hypgeom_bessel_i(iv, nu, z, prec);
        }

        acb_mul_onei(t, z);
        acb_hypgeom_bessel_j(jv, nu, t, prec);
        acb_pow(t, z, nu, prec);
        acb_mul(jv, jv, t, prec);
        acb_mul_onei(t, z);
        acb_pow(t, t, nu, prec);
        acb_div(jv, jv, t, prec);

        if (!acb_overlaps(iv, jv))
        {
            printf("FAIL: consistency with bessel_j\n\n");
            printf("nu = "); acb_printd(nu, 30); printf("\n\n");
            printf("z = "); acb_printd(z, 30); printf("\n\n");
            printf("iv = "); acb_printd(iv, 30); printf("\n\n");
            printf("jv = "); acb_printd(jv, 30); printf("\n\n");
            abort();
        }

        acb_clear(nu);
        acb_clear(z);
        acb_clear(jv);
        acb_clear(iv);
        acb_clear(t);
    }

    for (iter = 0; iter < 2000; iter++)
    {
        acb_t nu0, nu1, nu2, z, w0, w1, w2, t, u;
        slong prec0, prec1, prec2;

        acb_init(nu0);
        acb_init(nu1);
        acb_init(nu2);
        acb_init(z);
        acb_init(w0);
        acb_init(w1);
        acb_init(w2);
        acb_init(t);
        acb_init(u);

        prec0 = 2 + n_randint(state, 1000);
        prec1 = 2 + n_randint(state, 1000);
        prec2 = 2 + n_randint(state, 1000);

        acb_randtest_param(nu0, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(z, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(w0, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(w1, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));
        acb_randtest(w2, state, 1 + n_randint(state, 1000), 1 + n_randint(state, 100));

        acb_sub_ui(nu1, nu0, 1, prec0);
        acb_sub_ui(nu2, nu0, 2, prec0);

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_bessel_i_asymp(w0, nu0, z, prec0);
                break;
            case 1:
                acb_hypgeom_bessel_i_0f1(w0, nu0, z, prec0);
                break;
            default:
                acb_hypgeom_bessel_i(w0, nu0, z, prec0);
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_bessel_i_asymp(w1, nu0, z, prec1);
                break;
            case 1:
                acb_hypgeom_bessel_i_0f1(w1, nu0, z, prec1);
                break;
            default:
                acb_hypgeom_bessel_i(w1, nu0, z, prec1);
        }

        if (!acb_overlaps(w0, w1))
        {
            printf("FAIL: consistency\n\n");
            printf("nu = "); acb_printd(nu0, 30); printf("\n\n");
            printf("z = "); acb_printd(z, 30); printf("\n\n");
            printf("w0 = "); acb_printd(w0, 30); printf("\n\n");
            printf("w1 = "); acb_printd(w1, 30); printf("\n\n");
            abort();
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_bessel_i_asymp(w1, nu1, z, prec1);
                break;
            case 1:
                acb_hypgeom_bessel_i_0f1(w1, nu1, z, prec1);
                break;
            default:
                acb_hypgeom_bessel_i(w1, nu1, z, prec1);
        }

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_bessel_i_asymp(w2, nu2, z, prec2);
                break;
            case 1:
                acb_hypgeom_bessel_i_0f1(w2, nu2, z, prec2);
                break;
            default:
                acb_hypgeom_bessel_i(w2, nu2, z, prec2);
        }

        acb_mul(t, w1, nu1, prec0);
        acb_mul_2exp_si(t, t, 1);
        acb_submul(t, w2, z, prec0);
        acb_addmul(t, w0, z, prec0);

        if (!acb_contains_zero(t))
        {
            printf("FAIL: contiguous relation\n\n");
            printf("nu = "); acb_printd(nu0, 30); printf("\n\n");
            printf("z = ");  acb_printd(z, 30); printf("\n\n");
            printf("w0 = "); acb_printd(w0, 30); printf("\n\n");
            printf("w1 = "); acb_printd(w1, 30); printf("\n\n");
            printf("w2 = "); acb_printd(w2, 30); printf("\n\n");
            printf("t = "); acb_printd(t, 30); printf("\n\n");
            abort();
        }

        acb_neg(t, nu0);

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_bessel_i_asymp(w2, t, z, prec2);
                break;
            case 1:
                acb_hypgeom_bessel_i_0f1(w2, t, z, prec2);
                break;
            default:
                acb_hypgeom_bessel_i(w2, t, z, prec2);
        }

        acb_mul(w1, w1, w2, prec2);
        acb_neg(t, nu1);

        switch (n_randint(state, 3))
        {
            case 0:
                acb_hypgeom_bessel_i_asymp(w2, t, z, prec2);
                break;
            case 1:
                acb_hypgeom_bessel_i_0f1(w2, t, z, prec2);
                break;
            default:
                acb_hypgeom_bessel_i(w2, t, z, prec2);
        }

        acb_mul(w0, w0, w2, prec2);
        acb_sub(w0, w1, w0, prec2);

        acb_sin_pi(t, nu0, prec2);
        acb_const_pi(u, prec2);
        acb_mul(u, u, z, prec2);
        acb_div(t, t, u, prec2);
        acb_mul_2exp_si(t, t, 1);

        if (!acb_overlaps(w0, t))
        {
            printf("FAIL: wronskian\n\n");
            printf("nu = "); acb_printd(nu0, 30); printf("\n\n");
            printf("z = ");  acb_printd(z, 30); printf("\n\n");
            printf("w0 = "); acb_printd(w0, 30); printf("\n\n");
            printf("t = "); acb_printd(t, 30); printf("\n\n");
            abort();
        }

        acb_clear(nu0);
        acb_clear(nu1);
        acb_clear(nu2);
        acb_clear(z);
        acb_clear(w0);
        acb_clear(w1);
        acb_clear(w2);
        acb_clear(t);
        acb_clear(u);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

