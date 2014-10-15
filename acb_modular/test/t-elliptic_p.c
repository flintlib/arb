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

#include "acb_modular.h"

#define NUM_TESTS 6
#define EPS 1e-13

/* z, tau, p(z, tau) checked with Mathematica:
   N[{z, tau, WeierstrassP[z, WeierstrassInvariants[{1, tau}/2]]}, 20] */
const double testdata[NUM_TESTS][6] = {
    { 1.4142135623730950488, 1.7320508075688772935,
      2.2360679774997896964, 2.6457513110645905905,
      -3.3440423818409419071, 0.1143522489547216990 },
    { -3.0, -2.0, -7.0, 3.0, -3.3640051024505898616, 0.0 },
    { 0.1, 0.0, 0.6, 0.2, 99.322596953997138519, 2.258818666973309701 },
    { 0.0, 0.1, 0.6, 0.2, -99.322596953997138519, -2.258818666973309701 },
    { 0.5, 0.0, 0.333333333333333333, 20.0, 6.5797362673929057459, 0.0 },
    { 3.6666666666666667, 2014.0, -3.1415926535897932385, 0.1,
      -36.988356711748296440, -102.242185767588750178 }
};

static void
acb_set_dddd(acb_t z, double a, double ar, double b, double br)
{
    arf_set_d(arb_midref(acb_realref(z)), a);
    mag_set_d(arb_radref(acb_realref(z)), ar);
    arf_set_d(arb_midref(acb_imagref(z)), b);
    mag_set_d(arb_radref(acb_imagref(z)), br);
}

int main()
{
    long iter;
    flint_rand_t state;

    printf("elliptic_p....");
    fflush(stdout);

    flint_randinit(state);

    /* check test values */
    for (iter = 0; iter < 100; iter++)
    {
        long i;

        acb_t z, tau, p1, p2;

        acb_init(z);
        acb_init(tau);
        acb_init(p1);
        acb_init(p2);

        for (i = 0; i < NUM_TESTS; i++)
        {
            acb_set_dddd(z, testdata[i][0], 0.0, testdata[i][1], 0.0);
            acb_set_dddd(tau, testdata[i][2], 0.0, testdata[i][3], 0.0);
            acb_set_dddd(p2, testdata[i][4], EPS, testdata[i][5], EPS);

            acb_modular_elliptic_p(p1, z, tau, 2 + n_randint(state, 1000));

            if (!acb_overlaps(p1, p2))
            {
                printf("FAIL (test value)\n");
                printf("tau = "); acb_printd(tau, 15); printf("\n\n");
                printf("z = "); acb_printd(z, 15); printf("\n\n");
                printf("p1 = "); acb_printd(p1, 15); printf("\n\n");
                printf("p2 = "); acb_printd(p2, 15); printf("\n\n");
                abort();
            }
        }

        acb_clear(z);
        acb_clear(tau);
        acb_clear(p1);
        acb_clear(p2);
    }

    /* Test periods */
    for (iter = 0; iter < 2000; iter++)
    {
        acb_t tau, z1, z2, p1, p2;
        long m, n, e0, prec0, prec1, prec2;

        acb_init(tau);
        acb_init(z1);
        acb_init(z2);
        acb_init(p1);
        acb_init(p2);

        e0 = 1 + n_randint(state, 10);
        prec0 = 2 + n_randint(state, 1000);
        prec1 = 2 + n_randint(state, 1000);
        prec2 = 2 + n_randint(state, 1000);

        acb_randtest(tau, state, prec0, e0);
        if (arf_sgn(arb_midref(acb_imagref(tau))) < 0)
            acb_neg(tau, tau);

        acb_randtest(z1, state, prec0, e0);
        acb_randtest(p1, state, prec0, e0);
        acb_randtest(p2, state, prec0, e0);

        /* z2 = z1 + m + n*tau */
        m = n_randint(state, 10);
        n = n_randint(state, 10);
        acb_add_ui(z2, z1, m, prec0);
        acb_addmul_ui(z2, tau, n, prec0);

        acb_modular_elliptic_p(p1, z1, tau, prec1);
        acb_modular_elliptic_p(p2, z2, tau, prec2);

        if (!acb_overlaps(p1, p2))
        {
            printf("FAIL (overlap)\n");
            printf("tau = "); acb_printd(tau, 15); printf("\n\n");
            printf("z1 = "); acb_printd(z1, 15); printf("\n\n");
            printf("z2 = "); acb_printd(z2, 15); printf("\n\n");
            printf("p1 = "); acb_printd(p1, 15); printf("\n\n");
            printf("p2 = "); acb_printd(p2, 15); printf("\n\n");
            abort();
        }

        acb_modular_elliptic_p(z1, z1, tau, prec1);

        if (!acb_overlaps(z1, p1))
        {
            printf("FAIL (aliasing)\n");
            printf("tau = "); acb_printd(tau, 15); printf("\n\n");
            printf("z1 = "); acb_printd(z1, 15); printf("\n\n");
            printf("p1 = "); acb_printd(p1, 15); printf("\n\n");
            abort();
        }

        acb_clear(tau);
        acb_clear(z1);
        acb_clear(z2);
        acb_clear(p1);
        acb_clear(p2);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

