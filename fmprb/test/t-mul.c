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

#include "fmprb.h"

int
fmpr_close(const fmpr_t a, const fmpr_t b)
{
    fmpr_t t;
    int res1, res2;

    fmpr_init(t);

    fmpr_mul_ui(t, b, 257, FMPRB_RAD_PREC, FMPR_RND_UP);
    fmpr_mul_2exp_si(t, t, -8);
    res1 = fmpr_cmp(a, t) <= 0;

    fmpr_mul_ui(t, a, 257, FMPRB_RAD_PREC, FMPR_RND_UP);
    fmpr_mul_2exp_si(t, t, -8);
    res2 = fmpr_cmp(b, t) <= 0;

    fmpr_clear(t);

    return res1 && res2;
}

int main()
{
    long iter, iter2;
    flint_rand_t state;

    printf("mul....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        fmprb_t a, b, c;
        fmpq_t x, y, z;

        fmprb_init(a);
        fmprb_init(b);
        fmprb_init(c);

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        fmprb_randtest(a, state, 1 + n_randint(state, 200), 10);
        fmprb_randtest(b, state, 1 + n_randint(state, 200), 10);
        fmprb_randtest(c, state, 1 + n_randint(state, 200), 10);

        fmprb_get_rand_fmpq(x, state, a, 1 + n_randint(state, 200));
        fmprb_get_rand_fmpq(y, state, b, 1 + n_randint(state, 200));

        fmprb_mul(c, a, b, 2 + n_randint(state, 200));
        fmpq_mul(z, x, y);

        if (!fmprb_contains_fmpq(c, z))
        {
            printf("FAIL: containment\n\n");
            printf("a = "); fmprb_print(a); printf("\n\n");
            printf("x = "); fmpq_print(x); printf("\n\n");
            printf("b = "); fmprb_print(b); printf("\n\n");
            printf("y = "); fmpq_print(y); printf("\n\n");
            printf("c = "); fmprb_print(c); printf("\n\n");
            printf("z = "); fmpq_print(z); printf("\n\n");
            abort();
        }

        fmprb_clear(a);
        fmprb_clear(b);
        fmprb_clear(c);

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
    }

    /* aliasing of c and a */
    for (iter = 0; iter < 10000; iter++)
    {
        fmprb_t a, b;
        fmpq_t x, y, z;

        fmprb_init(a);
        fmprb_init(b);

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        fmprb_randtest(a, state, 1 + n_randint(state, 200), 10);
        fmprb_randtest(b, state, 1 + n_randint(state, 200), 10);

        fmprb_get_rand_fmpq(x, state, a, 1 + n_randint(state, 200));
        fmprb_get_rand_fmpq(y, state, b, 1 + n_randint(state, 200));

        fmprb_mul(a, a, b, 2 + n_randint(state, 200));
        fmpq_mul(z, x, y);

        if (!fmprb_contains_fmpq(a, z))
        {
            printf("FAIL: aliasing (c, a)\n\n");
            printf("a = "); fmprb_print(a); printf("\n\n");
            printf("x = "); fmpq_print(x); printf("\n\n");
            printf("b = "); fmprb_print(b); printf("\n\n");
            printf("y = "); fmpq_print(y); printf("\n\n");
            printf("z = "); fmpq_print(z); printf("\n\n");
            abort();
        }

        fmprb_clear(a);
        fmprb_clear(b);

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
    }

    /* aliasing of c and b */
    for (iter = 0; iter < 10000; iter++)
    {
        fmprb_t a, b;
        fmpq_t x, y, z;

        fmprb_init(a);
        fmprb_init(b);

        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        fmprb_randtest(a, state, 1 + n_randint(state, 200), 10);
        fmprb_randtest(b, state, 1 + n_randint(state, 200), 10);

        fmprb_get_rand_fmpq(x, state, a, 1 + n_randint(state, 200));
        fmprb_get_rand_fmpq(y, state, b, 1 + n_randint(state, 200));

        fmprb_mul(b, a, b, 2 + n_randint(state, 200));
        fmpq_mul(z, x, y);

        if (!fmprb_contains_fmpq(b, z))
        {
            printf("FAIL: aliasing (c, b)\n\n");
            printf("a = "); fmprb_print(a); printf("\n\n");
            printf("x = "); fmpq_print(x); printf("\n\n");
            printf("b = "); fmprb_print(b); printf("\n\n");
            printf("y = "); fmpq_print(y); printf("\n\n");
            printf("z = "); fmpq_print(z); printf("\n\n");
            abort();
        }

        fmprb_clear(a);
        fmprb_clear(b);

        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
    }

    /* main test */
    for (iter = 0; iter < 10000; iter++)
    {
        fmprb_t x, y, z, v;
        long prec;

        fmprb_init(x);
        fmprb_init(y);
        fmprb_init(z);
        fmprb_init(v);

        for (iter2 = 0; iter2 < 100; iter2++)
        {
            fmprb_randtest_special(x, state, n_randint(state,2) ? 2000 : 200, 200);
            fmprb_randtest_special(y, state, n_randint(state,2) ? 2000 : 200, 200);

            prec = 2 + n_randint(state, 2000);

            switch (n_randint(state, 5))
            {
            case 0:
                fmprb_mul(z, x, y, prec);
                fmprb_mul_naive(v, x, y, prec);

                if (!fmpr_equal(fmprb_midref(z), fmprb_midref(v))
                    || !fmpr_close(fmprb_radref(z), fmprb_radref(v)))
                {
                    printf("FAIL!\n");
                    printf("x = "); fmprb_print(x); printf("\n\n");
                    printf("y = "); fmprb_print(y); printf("\n\n");
                    printf("z = "); fmprb_print(z); printf("\n\n");
                    printf("v = "); fmprb_print(v); printf("\n\n");
                    abort();
                }
                break;

            case 1:
                fmprb_set(y, x);
                fmprb_mul(z, x, y, prec);
                fmprb_mul(v, x, x, prec);
                if (!fmprb_equal(z, v))
                {
                    printf("FAIL (aliasing 1)!\n");
                    printf("x = "); fmprb_print(x); printf("\n\n");
                    printf("z = "); fmprb_print(z); printf("\n\n");
                    printf("v = "); fmprb_print(v); printf("\n\n");
                    abort();
                }
                break;

            case 2:
                fmprb_mul(v, x, x, prec);
                fmprb_mul(x, x, x, prec);
                if (!fmprb_equal(v, x))
                {
                    printf("FAIL (aliasing 2)!\n");
                    printf("x = "); fmprb_print(x); printf("\n\n");
                    printf("z = "); fmprb_print(z); printf("\n\n");
                    printf("v = "); fmprb_print(v); printf("\n\n");
                    abort();
                }
                break;

            case 3:
                fmprb_mul(v, x, y, prec);
                fmprb_mul(x, x, y, prec);
                if (!fmprb_equal(x, v))
                {
                    printf("FAIL (aliasing 3)!\n");
                    printf("x = "); fmprb_print(x); printf("\n\n");
                    printf("y = "); fmprb_print(y); printf("\n\n");
                    printf("v = "); fmprb_print(v); printf("\n\n");
                    abort();
                }
                break;

            default:
                fmprb_mul(v, x, y, prec);
                fmprb_mul(x, y, x, prec);
                if (!fmprb_equal(x, v))
                {
                    printf("FAIL (aliasing 4)!\n");
                    printf("x = "); fmprb_print(x); printf("\n\n");
                    printf("y = "); fmprb_print(y); printf("\n\n");
                    printf("v = "); fmprb_print(v); printf("\n\n");
                    abort();
                }
                break;
            }
        }

        fmprb_clear(x);
        fmprb_clear(y);
        fmprb_clear(z);
        fmprb_clear(v);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
