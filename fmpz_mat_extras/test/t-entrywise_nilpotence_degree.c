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

#include "fmpz_mat_extras.h"
#include "perm.h"


int
_nilpotence_degree_is_superficially_ok_entrywise(const fmpz_mat_t A)
{
    slong n, i, j, d;

    if (!fmpz_mat_is_square(A))
        return 0;

    n = fmpz_mat_nrows(A);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            d = fmpz_get_si(fmpz_mat_entry(A, i, j));

            if (d < -1 || d > n)
                return 0;

            if (i == j && d == 0)
                return 0;

            if (i != j && d == 1)
                return 0;
        }
    }
    return 1;
}

/* todo: unfinished after here ... just copied from transitive closure ... */

/* permute rows and columns of a square matrix */
void
_fmpz_mat_permute(fmpz_mat_t B, const fmpz_mat_t A, const slong *perm)
{
    slong n, i, j;
    if (!fmpz_mat_is_square(A)) abort();
    if (A == B) abort();
    n = fmpz_mat_nrows(A);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            fmpz_set(
                    fmpz_mat_entry(B, perm[i], perm[j]),
                    fmpz_mat_entry(A, i, j));
        }
    }
}


/* this is not efficient */
void
_brute_force_entrywise_nilpotence_degree(fmpz_mat_t B, const fmpz_mat_t A)
{
    slong n, i, j, k;
    fmpz_mat_t S, curr, accum;

    n = fmpz_mat_nrows(A);
    fmpz_mat_init(S, n, n);
    fmpz_mat_init(curr, n, n);
    fmpz_mat_init(accum, n, n);
    fmpz_mat_entrywise_not_is_zero(S, A);
    fmpz_mat_one(curr);
    fmpz_mat_zero(accum);
    for (k = 0; k < 2*n+1; k++)
    {
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                if (!fmpz_is_zero(fmpz_mat_entry(curr, i, j)))
                {
                    fmpz_set_si(fmpz_mat_entry(accum, i, j), k+1);
                }
            }
        }
        fmpz_mat_mul(curr, curr, S);
    }
    fmpz_mat_clear(S);
    fmpz_mat_clear(curr);
    {
        slong i, j;
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                if (fmpz_cmp_si(fmpz_mat_entry(accum, i, j), n) > 0)
                {
                    fmpz_set_si(fmpz_mat_entry(B, i, j), -1);
                }
                else
                {
                    fmpz_set(fmpz_mat_entry(B, i, j),
                             fmpz_mat_entry(accum, i, j));

                }
            }
        }
    }
    fmpz_mat_clear(accum);
}


int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("entrywise_nilpotence_degree....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        slong m;
        fmpz_mat_t A, B, C, D;

        m = n_randint(state, 50);

        fmpz_mat_init(A, m, m);
        fmpz_mat_init(B, m, m);
        fmpz_mat_init(C, m, m);
        fmpz_mat_init(D, m, m);

        fmpz_mat_randtest(A, state, n_randint(state, 20) + 1);
        fmpz_mat_entrywise_not_is_zero(A, A);

        fmpz_mat_entrywise_nilpotence_degree(B, A);

        if (!_nilpotence_degree_is_superficially_ok_entrywise(B))
        {
            flint_printf("FAIL (entrywise)\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(B); flint_printf("\n\n");
            abort();
        }

        /* aliasing */
        {
            fmpz_mat_set(C, A);
            fmpz_mat_entrywise_nilpotence_degree(C, C);
            if (!fmpz_mat_equal(B, C))
            {
                flint_printf("FAIL (aliasing)\n");
                fmpz_mat_print_pretty(A); flint_printf("\n\n");
                fmpz_mat_print_pretty(B); flint_printf("\n\n");
                fmpz_mat_print_pretty(C); flint_printf("\n\n");
                abort();
            }
        }

        /* reduction to transitive closure */
        {
            slong i, j;
            fmpz_mat_t U, V;

            fmpz_mat_init(U, m, m);
            fmpz_mat_transitive_closure(U, A);
            fmpz_mat_entrywise_not_is_zero(U, U);

            fmpz_mat_init(V, m, m);
            for (i = 0; i < m; i++)
            {
                for (j = 0; j < m; j++)
                {
                    if (fmpz_is_zero(fmpz_mat_entry(B, i, j)) ||
                        fmpz_is_one(fmpz_mat_entry(B, i, j)))
                    {
                        fmpz_zero(fmpz_mat_entry(V, i, j));
                    }
                    else
                    {
                        fmpz_one(fmpz_mat_entry(V, i, j));
                    }
                }
            }

            if (!fmpz_mat_equal(U, V))
            {
                flint_printf("FAIL (reduction to transitive closure)\n");
                fmpz_mat_print_pretty(A); flint_printf("\n\n");
                fmpz_mat_print_pretty(B); flint_printf("\n\n");
                fmpz_mat_print_pretty(U); flint_printf("\n\n");
                fmpz_mat_print_pretty(V); flint_printf("\n\n");
                abort();
            }

            fmpz_mat_clear(U);
            fmpz_mat_clear(V);
        }

        /* test commutativity of permutation with entrywise nilpotence */
        {
            slong *perm;
            perm = flint_malloc(m * sizeof(slong));
            _perm_randtest(perm, m, state);

            /* C is the entrywise nilpotence of the permutation of A */
            fmpz_mat_randtest(C, state, n_randint(state, 20) + 1);
            _fmpz_mat_permute(C, A, perm);
            fmpz_mat_entrywise_nilpotence_degree(C, C);

            /* D is the permutation of the entrywise nilpotence of A */
            fmpz_mat_randtest(D, state, n_randint(state, 20) + 1);
            _fmpz_mat_permute(D, B, perm);

            if (!fmpz_mat_equal(C, D))
            {
                flint_printf("FAIL (commutativity with permutation)\n");
                fmpz_mat_print_pretty(A); flint_printf("\n\n");
                fmpz_mat_print_pretty(B); flint_printf("\n\n");
                fmpz_mat_print_pretty(C); flint_printf("\n\n");
                fmpz_mat_print_pretty(D); flint_printf("\n\n");
                abort();
            }
            flint_free(perm);
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
        fmpz_mat_clear(D);
    }

    /* use brute force to check small examples */
    for (iter = 0; iter < 1000; iter++)
    {
        slong m;
        fmpz_mat_t A, B, C;

        m = n_randint(state, 10);

        fmpz_mat_init(A, m, m);
        fmpz_mat_init(B, m, m);
        fmpz_mat_init(C, m, m);

        fmpz_mat_randtest(A, state, n_randint(state, 20) + 1);
        fmpz_mat_entrywise_not_is_zero(A, A);

        fmpz_mat_entrywise_nilpotence_degree(B, A);
        _brute_force_entrywise_nilpotence_degree(C, A);

        if (!fmpz_mat_equal(B, C))
        {
            flint_printf("FAIL (brute force)\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(B); flint_printf("\n\n");
            fmpz_mat_print_pretty(C); flint_printf("\n\n");
            abort();
        }

        fmpz_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_mat_clear(C);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
