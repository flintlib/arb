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
#include "arith.h"

static void
acb_modular_epsilon_arg_naive(fmpq_t arg, const psl2z_t g)
{
#define a (&g->a)
#define b (&g->b)
#define c (&g->c)
#define d (&g->d)

    if (fmpz_is_zero(c))
    {
        fmpz_set(fmpq_numref(arg), b);
        fmpz_set_ui(fmpq_denref(arg), 12);
        fmpq_canonicalise(arg);
    }
    else
    {
        fmpq_t t;
        fmpq_init(t);

        fmpz_add(fmpq_numref(arg), a, d);
        fmpz_submul_ui(fmpq_numref(arg), c, 3);
        fmpz_mul_ui(fmpq_denref(arg), c, 12);
        fmpq_canonicalise(arg);

        arith_dedekind_sum(t, d, c);
        fmpq_sub(arg, arg, t);

        fmpq_clear(t);
    }
#undef a
#undef b
#undef c
#undef d
}

int main()
{
    long iter;
    flint_rand_t state;

    printf("eta_epsilon_arg....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        psl2z_t g;
        fmpq_t x, y;
        acb_t a, b;

        psl2z_init(g);
        fmpq_init(x);
        fmpq_init(y);
        acb_init(a);
        acb_init(b);

        psl2z_randtest(g, state, n_randint(state, 200));
        fmpq_randtest(x, state, 100);
        fmpq_randtest(y, state, 100);

        fmpq_set_si(x, acb_modular_epsilon_arg(g), 12);
        acb_modular_epsilon_arg_naive(y, g);

        arb_sin_cos_pi_fmpq(acb_imagref(a), acb_realref(a), x, 200);
        arb_sin_cos_pi_fmpq(acb_imagref(b), acb_realref(b), x, 200);

        if (!acb_overlaps(a, b))
        {
            printf("FAIL\n");
            printf("g = "); psl2z_print(g); printf("\n");
            printf("x = "); fmpq_print(x); printf("\n");
            printf("y = "); fmpq_print(y); printf("\n");
            abort();
        }

        psl2z_clear(g);
        fmpq_clear(x);
        fmpq_clear(y);
        acb_clear(a);
        acb_clear(b);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

