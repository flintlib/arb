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

#include "arb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("exp_sum_bs_powtab....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        long N;
        fmpz_t x, T, Q;
        fmpq_t S, V;
        mp_bitcnt_t Qexp, r;

        fmpz_init(x);
        fmpz_init(T);
        fmpz_init(Q);
        fmpq_init(S);
        fmpq_init(V);

        N = 1 + n_randint(state, 300);
        r = n_randint(state, 10);
        fmpz_randtest(x, state, 80);

        _arb_exp_sum_bs_simple(T, Q, &Qexp, x, r, N);
        fmpq_set_fmpz_frac(S, T, Q);
        fmpq_div_2exp(S, S, Qexp);

        _arb_exp_sum_bs_powtab(T, Q, &Qexp, x, r, N);
        fmpq_set_fmpz_frac(V, T, Q);
        fmpq_div_2exp(V, V, Qexp);

        if (!fmpq_equal(S, V))
        {
            printf("FAIL\n\n");
            printf("N = %ld\n\n", N);
            printf("r = %lu\n\n", r);
            printf("x = "); fmpz_print(x); printf("\n\n");
            printf("T = "); fmpz_print(T); printf("\n\n");
            printf("Q = "); fmpz_print(T); printf("\n\n");
            printf("V = "); fmpq_print(V); printf("\n\n");
            printf("S = "); fmpq_print(S); printf("\n\n");
            abort();
        }

        fmpz_clear(x);
        fmpz_clear(T);
        fmpz_clear(Q);
        fmpq_clear(S);
        fmpq_clear(V);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

