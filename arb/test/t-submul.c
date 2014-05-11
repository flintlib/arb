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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "arb.h"

int
mag_close(const mag_t am, const mag_t bm)
{
    arf_t t, a, b;
    int res1, res2;

    arf_init(t);
    arf_init(a);
    arf_init(b);

    arf_set_mag(a, am);
    arf_set_mag(b, bm);

    arf_mul_ui(t, b, 257, MAG_BITS, ARF_RND_UP);
    arf_mul_2exp_si(t, t, -8);
    res1 = arf_cmp(a, t) <= 0;

    arf_mul_ui(t, a, 257, MAG_BITS, ARF_RND_UP);
    arf_mul_2exp_si(t, t, -8);
    res2 = arf_cmp(b, t) <= 0;

    arf_clear(t);
    arf_clear(a);
    arf_clear(b);

    return res1 && res2;
}

void
arb_submul_naive(arb_t z, const arb_t x, const arb_t y, long prec)
{
    arb_t t;
    arb_init(t);
    arb_mul(t, x, y, ARF_PREC_EXACT);
    arb_sub(z, z, t, prec);
    arb_clear(t);
}

int main()
{
    long iter, iter2;
    flint_rand_t state;

    printf("submul....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        arb_t a, b, c;
        fmpq_t x, y, z;

        arb_init(a);
        arb_init(b);
        arb_init(c);

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        arb_randtest(a, state, 1 + n_randint(state, 200), 10);
        arb_randtest(b, state, 1 + n_randint(state, 200), 10);
        arb_randtest(c, state, 1 + n_randint(state, 200), 10);

        arb_get_rand_fmpq(x, state, a, 1 + n_randint(state, 200));
        arb_get_rand_fmpq(y, state, b, 1 + n_randint(state, 200));
        arb_get_rand_fmpq(z, state, c, 1 + n_randint(state, 200));

        arb_submul(c, a, b, 2 + n_randint(state, 200));
        fmpq_submul(z, x, y);

        if (!arb_contains_fmpq(c, z))
        {
            printf("FAIL: containment\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("x = "); fmpq_print(x); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            printf("y = "); fmpq_print(y); printf("\n\n");
            printf("c = "); arb_print(c); printf("\n\n");
            printf("z = "); fmpq_print(z); printf("\n\n");
            abort();
        }

        arb_clear(a);
        arb_clear(b);
        arb_clear(c);

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
    }

    /* aliasing of c and a */
    for (iter = 0; iter < 10000; iter++)
    {
        arb_t a, b;
        fmpq_t x, y, z;

        arb_init(a);
        arb_init(b);

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        arb_randtest(a, state, 1 + n_randint(state, 200), 10);
        arb_randtest(b, state, 1 + n_randint(state, 200), 10);

        arb_get_rand_fmpq(x, state, a, 1 + n_randint(state, 200));
        arb_get_rand_fmpq(y, state, b, 1 + n_randint(state, 200));
        fmpq_set(z, x);

        arb_submul(a, a, b, 2 + n_randint(state, 200));
        fmpq_submul(z, x, y);

        if (!arb_contains_fmpq(a, z))
        {
            printf("FAIL: aliasing (c, a)\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("x = "); fmpq_print(x); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            printf("y = "); fmpq_print(y); printf("\n\n");
            printf("z = "); fmpq_print(z); printf("\n\n");
            abort();
        }

        arb_clear(a);
        arb_clear(b);

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
    }

    /* aliasing of c and b */
    for (iter = 0; iter < 10000; iter++)
    {
        arb_t a, b;
        fmpq_t x, y, z;

        arb_init(a);
        arb_init(b);

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        arb_randtest(a, state, 1 + n_randint(state, 200), 10);
        arb_randtest(b, state, 1 + n_randint(state, 200), 10);

        arb_get_rand_fmpq(x, state, a, 1 + n_randint(state, 200));
        arb_get_rand_fmpq(y, state, b, 1 + n_randint(state, 200));
        fmpq_set(z, y);

        arb_submul(b, a, b, 2 + n_randint(state, 200));
        fmpq_submul(z, x, y);

        if (!arb_contains_fmpq(b, z))
        {
            printf("FAIL: aliasing (c, b)\n\n");
            printf("a = "); arb_print(a); printf("\n\n");
            printf("x = "); fmpq_print(x); printf("\n\n");
            printf("b = "); arb_print(b); printf("\n\n");
            printf("y = "); fmpq_print(y); printf("\n\n");
            printf("z = "); fmpq_print(z); printf("\n\n");
            abort();
        }

        arb_clear(a);
        arb_clear(b);

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
    }

    /* main test */
    for (iter = 0; iter < 10000; iter++)
    {
        arb_t x, y, z, v;
        long prec;

        arb_init(x);
        arb_init(y);
        arb_init(z);
        arb_init(v);

        for (iter2 = 0; iter2 < 100; iter2++)
        {
            arb_randtest_special(x, state, n_randint(state,2) ? 2000 : 200, 200);
            arb_randtest_special(y, state, n_randint(state,2) ? 2000 : 200, 200);
            arb_randtest_special(z, state, n_randint(state,2) ? 2000 : 200, 200);

            prec = 2 + n_randint(state, 2000);

            switch (n_randint(state, 5))
            {
            case 0:
                arb_set(v, z);
                arb_submul(z, x, y, prec);
                arb_submul_naive(v, x, y, prec);

                if (!arf_equal(arb_midref(z), arb_midref(v))
                    || !mag_close(arb_radref(z), arb_radref(v)))
                {
                    printf("FAIL!\n");
                    printf("x = "); arb_print(x); printf("\n\n");
                    printf("y = "); arb_print(y); printf("\n\n");
                    printf("z = "); arb_print(z); printf("\n\n");
                    printf("v = "); arb_print(v); printf("\n\n");
                    abort();
                }
                break;

            case 1:
                arb_set(y, x);
                arb_set(z, v);
                arb_submul(z, x, y, prec);
                arb_submul(v, x, x, prec);

                if (!arf_equal(arb_midref(z), arb_midref(v))
                    || !mag_close(arb_radref(z), arb_radref(v)))
                {
                    printf("FAIL (aliasing 1)!\n");
                    printf("x = "); arb_print(x); printf("\n\n");
                    printf("z = "); arb_print(z); printf("\n\n");
                    printf("v = "); arb_print(v); printf("\n\n");
                    abort();
                }
                break;

            case 2:
                arb_set(v, x);
                arb_submul(v, x, x, prec);
                arb_submul(x, x, x, prec);

                if (!arf_equal(arb_midref(x), arb_midref(v))
                    || !mag_close(arb_radref(x), arb_radref(v)))
                {
                    printf("FAIL (aliasing 2)!\n");
                    printf("x = "); arb_print(x); printf("\n\n");
                    printf("v = "); arb_print(v); printf("\n\n");
                    abort();
                }
                break;

            case 3:
                arb_set(v, x);
                arb_submul(v, x, y, prec);
                arb_submul(x, x, y, prec);

                if (!arf_equal(arb_midref(x), arb_midref(v))
                    || !mag_close(arb_radref(x), arb_radref(v)))
                {
                    printf("FAIL (aliasing 3)!\n");
                    printf("x = "); arb_print(x); printf("\n\n");
                    printf("y = "); arb_print(y); printf("\n\n");
                    printf("v = "); arb_print(v); printf("\n\n");
                    abort();
                }
                break;

            default:
                arb_set(v, x);
                arb_submul(v, x, y, prec);
                arb_submul(x, y, x, prec);

                if (!arf_equal(arb_midref(x), arb_midref(v))
                    || !mag_close(arb_radref(x), arb_radref(v)))
                {
                    printf("FAIL (aliasing 4)!\n");
                    printf("x = "); arb_print(x); printf("\n\n");
                    printf("y = "); arb_print(y); printf("\n\n");
                    printf("v = "); arb_print(v); printf("\n\n");
                    abort();
                }
                break;
            }
        }

        arb_clear(x);
        arb_clear(y);
        arb_clear(z);
        arb_clear(v);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
