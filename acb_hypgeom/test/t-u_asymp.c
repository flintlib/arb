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

void
acb_randtest_maybe_half_int(acb_t x, flint_rand_t state, long prec, long size)
{
    if (n_randint(state, 8) == 0)
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_randtest(t, state, 1 + n_randint(state, prec));
        arb_set_fmpz(acb_realref(x), t);
        arb_zero(acb_imagref(x));
        acb_mul_2exp_si(x, x, -1);
        fmpz_clear(t);
    }
    else
    {
        acb_randtest(x, state, prec, size);
    }
}

int main()
{
    long iter;
    flint_rand_t state;

    printf("u_asymp....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        acb_t a, b, a2, b2, z, U1, U2, t, u, M1, M2, am;
        acb_struct bm[2];
        ulong n1, n2;
        long prec0, prec1, prec2;

        acb_init(a);
        acb_init(b);
        acb_init(a2);
        acb_init(b2);
        acb_init(z);
        acb_init(U1);
        acb_init(U2);
        acb_init(t);
        acb_init(u);
        acb_init(M1);
        acb_init(M2);
        acb_init(am);
        acb_init(bm);
        acb_init(bm + 1);

        if (n_randint(state, 4) == 0)
        {
            n1 = n_randint(state, 200);
            n2 = n_randint(state, 200);
            prec0 = 2 + n_randint(state, 1000);
            prec1 = 2 + n_randint(state, 1000);
            prec2 = 2 + n_randint(state, 1000);
        }
        else
        {
            n1 = n_randint(state, 40);
            n2 = n_randint(state, 40);
            prec0 = 2 + n_randint(state, 300);
            prec1 = 2 + n_randint(state, 300);
            prec2 = 2 + n_randint(state, 300);
        }

        acb_randtest_maybe_half_int(a, state, prec0, 1 + n_randint(state, 20));

        if (n_randint(state, 4) == 0)
            acb_add_ui(b, a, n_randint(state, 10), prec0);
        else
            acb_randtest_maybe_half_int(b, state, prec0, 1 + n_randint(state, 20));

        acb_randtest(z, state, prec0, 1 + n_randint(state, 20));

        /* Test Kummer's transformation */
        acb_sub(a2, a, b, prec0);
        acb_add_ui(a2, a2, 1, prec0);
        acb_sub_ui(b2, b, 2, prec0);
        acb_neg(b2, b2);

        acb_hypgeom_u_asymp(U1, a, b, z, n1, prec1);
        acb_hypgeom_u_asymp(U2, a2, b2, z, n2, prec2);

        if (!acb_overlaps(U1, U2))
        {
            printf("FAIL (Kummer transformation)\n");
            printf("iter = %ld\n", iter);
            printf("a = "); acb_printd(a, 50); printf("\n");
            printf("b = "); acb_printd(b, 50); printf("\n");
            printf("z = "); acb_printd(z, 50); printf("\n");
            printf("n1 = %ld, n2 = %ld, prec1 = %ld, prec2 = %ld\n", n1, n2, prec1, prec2);
            printf("U1 = "); acb_printd(U1, 100); printf("\n");
            printf("U2 = "); acb_printd(U2, 100); printf("\n");
            abort();
        }

        /* Check contiguous relation
            (b-a-1)U(a,b-1,z) + z U(a,b+1,z) + (1-b-z) U(a,b,z) = 0 */
        acb_one(t);
        acb_sub(t, t, b, prec0);
        acb_sub(t, t, z, prec0);
        acb_mul(u, U1, t, prec1);

        acb_add_ui(b2, b, 1, prec0);
        acb_hypgeom_u_asymp(U2, a, b2, z, n2, prec2);
        acb_addmul(u, U2, z, prec1);

        acb_sub_ui(b2, b, 1, prec0);
        acb_hypgeom_u_asymp(U2, a, b2, z, n2, prec2);
        acb_sub(t, b, a, prec0);
        acb_sub_ui(t, t, 1, prec0);
        acb_mul(t, t, U2, prec0);
        acb_add(t, t, u, prec0);

        if (!acb_contains_zero(t))
        {
            printf("FAIL (contiguous relation)\n");
            printf("iter = %ld\n", iter);
            printf("a = "); acb_printd(a, 50); printf("\n");
            printf("b = "); acb_printd(b, 50); printf("\n");
            printf("z = "); acb_printd(z, 50); printf("\n");
            printf("n1 = %ld, n2 = %ld, prec1 = %ld, prec2 = %ld\n", n1, n2, prec1, prec2);
            printf("U1 = "); acb_printd(U1, 100); printf("\n");
            printf("t = "); acb_printd(t, 100); printf("\n");
            abort();
        }

        /* Test U(a,b,z) = gamma(1-b)/gamma(a-b+1) M(a,b,z)
                         + gamma(b-1)/gamma(a) z^(1-b) M(a-b+1,2-b,z) */
        acb_set(am, a);
        acb_set(bm, b);
        acb_one(bm + 1);
        acb_hypgeom_pfq_direct(M1, am, 1, bm, 2, z, n2, prec2);

        acb_sub(am, a, b, prec2);
        acb_add_ui(am, am, 1, prec2);
        acb_sub_ui(bm, b, 2, prec2);
        acb_neg(bm, bm);
        acb_one(bm + 1);
        acb_hypgeom_pfq_direct(M2, am, 1, bm, 2, z, n2, prec2);

        acb_sub(am, a, b, prec2);
        acb_add_ui(am, am, 1, prec2);
        acb_rgamma(am, am, prec2);
        acb_mul(M1, M1, am, prec2);
        acb_sub_ui(am, b, 1, prec2);
        acb_neg(am, am);
        acb_gamma(am, am, prec2);
        acb_mul(M1, M1, am, prec2);

        acb_rgamma(am, a, prec2);
        acb_mul(M2, M2, am, prec2);
        acb_sub_ui(am, b, 1, prec2);
        acb_gamma(am, am, prec2);
        acb_mul(M2, M2, am, prec2);

        acb_sub_ui(am, b, 1, prec2);
        acb_neg(am, am);
        acb_pow(am, z, am, prec2);
        acb_mul(M2, M2, am, prec2);

        acb_add(U2, M1, M2, prec2);

        acb_pow(am, z, a, prec2);
        acb_mul(U2, U2, am, prec2);

        if (!acb_overlaps(U1, U2))
        {
            printf("FAIL (U in terms of M)\n");
            printf("iter = %ld\n", iter);
            printf("a = "); acb_printd(a, 50); printf("\n");
            printf("b = "); acb_printd(b, 50); printf("\n");
            printf("z = "); acb_printd(z, 50); printf("\n");
            printf("n1 = %ld, n2 = %ld, prec1 = %ld, prec2 = %ld\n", n1, n2, prec1, prec2);
            printf("U1 = "); acb_printd(U1, 100); printf("\n");
            printf("U2 = "); acb_printd(U2, 100); printf("\n");
            abort();
        }

        /* Test special value: b = a+1 */
        acb_add_ui(b, a, 1, prec0);
        acb_hypgeom_u_asymp(U1, a, b, z, n1, prec1);
        acb_one(U2);

        if (!acb_overlaps(U1, U2))
        {
            printf("FAIL (special value)\n");
            printf("iter = %ld\n", iter);
            printf("a = "); acb_printd(a, 50); printf("\n");
            printf("b = "); acb_printd(b, 50); printf("\n");
            printf("z = "); acb_printd(z, 50); printf("\n");
            printf("n1 = %ld, n2 = %ld, prec1 = %ld, prec2 = %ld\n", n1, n2, prec1, prec2);
            printf("U1 = "); acb_printd(U1, 100); printf("\n");
            printf("U2 = "); acb_printd(U2, 100); printf("\n");
            abort();
        }

        acb_clear(a);
        acb_clear(b);
        acb_clear(a2);
        acb_clear(b2);
        acb_clear(z);
        acb_clear(U1);
        acb_clear(U2);
        acb_clear(t);
        acb_clear(u);
        acb_clear(M1);
        acb_clear(M2);
        acb_clear(am);
        acb_clear(bm);
        acb_clear(bm + 1);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

