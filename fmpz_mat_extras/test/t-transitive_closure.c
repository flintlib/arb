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

/* transitive closure can only turn zeros into ones */
int
_transitive_closure_is_ok_entrywise(const fmpz_mat_t X, const fmpz_mat_t Y)
{
    slong i, j;
    if (fmpz_mat_nrows(X) != fmpz_mat_nrows(Y) ||
        fmpz_mat_ncols(X) != fmpz_mat_ncols(Y))
    {
        return 0;
    }
    for (i = 0; i < fmpz_mat_nrows(X); i++)
    {
        for (j = 0; j < fmpz_mat_ncols(X); j++)
        {
            if (!fmpz_equal(
                        fmpz_mat_entry(X, i, j),
                        fmpz_mat_entry(Y, i, j)))
            {
                if (!fmpz_is_zero(fmpz_mat_entry(X, i, j)))
                {
                    return 0;
                }
                if (!fmpz_is_one(fmpz_mat_entry(Y, i, j)))
                {
                    return 0;
                }
            }
        }
    }
    return 1;
}

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
_brute_force_transitive_closure(fmpz_mat_t B, const fmpz_mat_t A)
{
    slong n, k;
    fmpz_mat_t S, curr, accum;

    n = fmpz_mat_nrows(A);
    fmpz_mat_init(S, n, n);
    fmpz_mat_init(curr, n, n);
    fmpz_mat_init(accum, n, n);
    fmpz_mat_entrywise_not_is_zero(S, A);
    fmpz_mat_one(curr);
    fmpz_mat_zero(accum);
    for (k = 0; k < n; k++)
    {
        fmpz_mat_mul(curr, curr, S);
        fmpz_mat_add(accum, accum, curr);
    }
    fmpz_mat_clear(S);
    fmpz_mat_clear(curr);
    {
        slong i, j;
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                if (fmpz_is_zero(fmpz_mat_entry(A, i, j)) &&
                    !fmpz_is_zero(fmpz_mat_entry(accum, i, j)))
                {
                    fmpz_one(fmpz_mat_entry(B, i, j));
                }
                else
                {
                    fmpz_set(fmpz_mat_entry(B, i, j),
                             fmpz_mat_entry(A, i, j));
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

    flint_printf("transitive_closure....");
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
        fmpz_mat_randtest(B, state, n_randint(state, 20) + 1);

        fmpz_mat_transitive_closure(B, A);

        /* test local properties of the closure */
        if (!_transitive_closure_is_ok_entrywise(A, B))
        {
            flint_printf("FAIL (entrywise)\n");
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(B); flint_printf("\n\n");
            abort();
        }

        /* test aliasing */
        {
            fmpz_mat_set(C, A);
            fmpz_mat_transitive_closure(C, C);
            if (!fmpz_mat_equal(B, C))
            {
                flint_printf("FAIL (aliasing)\n");
                fmpz_mat_print_pretty(A); flint_printf("\n\n");
                fmpz_mat_print_pretty(B); flint_printf("\n\n");
                fmpz_mat_print_pretty(C); flint_printf("\n\n");
                abort();
            }
        }

        /* test commutativity of permutation with transitive closure */
        {
            slong *perm;
            perm = flint_malloc(m * sizeof(slong));
            _perm_randtest(perm, m, state);

            /* C is the transitive closure of the permutation of A */
            fmpz_mat_randtest(C, state, n_randint(state, 20) + 1);
            _fmpz_mat_permute(C, A, perm);
            fmpz_mat_transitive_closure(C, C);

            /* D is the permutation of the transitive closure of A */
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

    /* check transitive closure using brute force with smallish matrices */
    for (iter = 0; iter < 1000; iter++)
    {
        slong m;
        fmpz_mat_t A, B, C;

        m = n_randint(state, 10);

        fmpz_mat_init(A, m, m);
        fmpz_mat_init(B, m, m);
        fmpz_mat_init(C, m, m);

        fmpz_mat_randtest(A, state, n_randint(state, 20) + 1);
        fmpz_mat_randtest(B, state, n_randint(state, 20) + 1);
        fmpz_mat_randtest(C, state, n_randint(state, 20) + 1);

        fmpz_mat_transitive_closure(B, A);
        _brute_force_transitive_closure(C, A);

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
