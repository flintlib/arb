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

#include "flint/fmpq_mat.h"
#include "arb_mat.h"

/*
 *  1 -1
 * -1  2 -1
 *    -1  2 -1
 *        .  .  .
 *          -1  2 -1
 *             -1  1
 */
static void
_dct2_A(arb_mat_t A)
{
    slong k, n;
    n = arb_mat_nrows(A);
    arb_mat_zero(A);
    for (k = 0; k < n-1; k++)
    {
        arb_set_si(arb_mat_entry(A, k+1, k), -1);
        arb_set_si(arb_mat_entry(A, k, k+1), -1);
    }
    for (k = 1; k < n-1; k++)
    {
        arb_set_si(arb_mat_entry(A, k, k), 2);
    }
    arb_set_si(arb_mat_entry(A, 0, 0), 1);
    arb_set_si(arb_mat_entry(A, n-1, n-1), 1);
}

void
_dct2_P(arb_mat_t P, slong prec)
{
    slong j, k, n;
    arb_t x, r;

    n = arb_mat_nrows(P);
    arb_init(x);
    arb_init(r);

    /* set the first column to 1/sqrt(n) */
    arb_set_si(x, n);
    arb_rsqrt(x, x, prec);
    for (j = 0; j < n; j++)
    {
        arb_set(arb_mat_entry(P, j, 0), x);
    }

    /* set remaining columns to sqrt(2/n)*cos[(j+1/2)*k*pi/n] */
    arb_set_si(r, n);
    arb_mul_2exp_si(r, r, -1);
    arb_rsqrt(r, r, prec);
    for (k = 1; k < n; k++)
    {
        for (j = 0; j < n; j++)
        {
            arb_set_d(x, (2*j + 1) * k);
            arb_div_si(x, x, 2*n, prec);
            arb_cos_pi(x, x, prec);
            arb_mul(arb_mat_entry(P, j, k), x, r, prec);
        }
    }
    arb_clear(x);
    arb_clear(r);
}

void
_dct2_D(arb_mat_t D, slong prec)
{
    slong k, n;
    arb_t x;

    n = arb_mat_nrows(D);
    arb_init(x);
    for (k = 0; k < n; k++)
    {
        arb_set_si(x, k);
        arb_div_si(x, x, n, prec);
        arb_cos_pi(x, x, prec);
        arb_sub_si(x, x, 1, prec);
        arb_neg(x, x);
        arb_mul_2exp_si(x, x, 1);
        arb_set(arb_mat_entry(D, k, 0), x);
    }
    arb_clear(x);
}

static void
_test_dct2(flint_rand_t state)
{
    slong n, p1, p2, logscale;

    /* n = n_randint(state, 8) + 2; */
    n = 5;
    p1 = 2 + n_randint(state, 202);
    p2 = 2 + n_randint(state, 202);

    /* this works for me */
    /* logscale = n_randint(state, 4); */

    /* but this one, not so much... */
    logscale = n_randint(state, 10);

    if (n_randint(state, 2))
        logscale = -logscale;

    {
        arb_mat_t A, D, P, Dt, Pt;

        arb_mat_init(A, n, n);
        arb_mat_init(D, n, 1);
        arb_mat_init(P, n, n);
        arb_mat_init(Dt, n, 1);
        arb_mat_init(Pt, n, n);

        /* compute the diagonalization of the scaled DCT-II matrix */
        _dct2_A(A);
        arb_mat_scalar_mul_2exp_si(A, A, logscale);
        arb_mat_symmetric_diagonalization(D, P, A, p1);

        /* compute the true eigenvectors and the true scaled eigenvalues */
        _dct2_D(Dt, p2);
        _dct2_P(Pt, p2);
        arb_mat_scalar_mul_2exp_si(Dt, Dt, logscale);

        if (!arb_mat_overlaps(D, Dt))
        {
            flint_printf("FAIL (DCT-II eigenvalues)\n");
            flint_printf("logscale=%wd p1=%wd p2=%wd\n", logscale, p1, p2);
            flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n");
            flint_printf("eigenvalues computed using diagonalization = \n");
            arb_mat_printd(D, 15); flint_printf("\n");
            flint_printf("true eigenvalues = \n");
            arb_mat_printd(Dt, 15); flint_printf("\n");
            abort();
        }

        if (!arb_mat_overlaps(P, Pt))
        {
            flint_printf("FAIL (DCT-II eigenvectors)\n");
            flint_printf("logscale=%wd p1=%wd p2=%wd\n", logscale, p1, p2);
            flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n");
            flint_printf("eigenvalues computed using diagonalization = \n");
            arb_mat_printd(D, 15); flint_printf("\n");
            flint_printf("true eigenvalues = \n");
            arb_mat_printd(Dt, 15); flint_printf("\n");
            flint_printf("eigenvectors computed using diagonalization = \n");
            arb_mat_printd(P, 15); flint_printf("\n");
            flint_printf("true eigenvectors = \n");
            arb_mat_printd(Pt, 15); flint_printf("\n");
            abort();
        }

        arb_mat_clear(A);
        arb_mat_clear(D);
        arb_mat_clear(P);
        arb_mat_clear(Dt);
        arb_mat_clear(Pt);
    }
}

static int
_fmpq_mat_is_symmetric(const fmpq_mat_t Q)
{
    slong i, j;
    if (!fmpq_mat_is_square(Q))
        return 0;
    for (i = 0; i < fmpq_mat_nrows(Q); i++)
        for (j = 0; j < i; j++)
            if (!fmpq_equal(fmpq_mat_entry(Q, i, j), fmpq_mat_entry(Q, j, i)))
                return 0;
    return 1;
}

static int
_fmpq_mat_is_one(const fmpq_mat_t Q)
{
    slong i, j, n;
    if (!fmpq_mat_is_square(Q))
        return 0;
    n = fmpq_mat_nrows(Q);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            fmpq *q = fmpq_mat_entry(Q, i, j);
            if (i == j && !fmpq_is_one(q))
                return 0;
            if (i != j && !fmpq_is_zero(q))
                return 0;
        }
    }
    return 1;
}

static void
_fmpq_mat_randtest_skew_symmetric(
        fmpq_mat_t Q, flint_rand_t state, mp_bitcnt_t bits)
{
    slong i, j, n;
    if (!fmpq_mat_is_square(Q)) abort(); /* assert */
    fmpq_mat_randtest(Q, state, bits);
    n = fmpq_mat_nrows(Q);
    for (i = 0; i < n; i++)
        fmpq_zero(fmpq_mat_entry(Q, i, i));
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < i; j++)
        {
            fmpq_neg(fmpq_mat_entry(Q, i, j),
                     fmpq_mat_entry(Q, j, i));
        }
    }
}

void
_fmpq_mat_randtest_orthogonal(
        fmpq_mat_t Q, flint_rand_t state, mp_bitcnt_t bits)
{
    /*
     * S -> (S - I)^-1 (S + I)
     * The Generation of All Rational Orthogonal Matrices
     * Hans Liebeck and Anthony Osborne
     * The American Mathematical Monthly,
     * Vol 98, No 2, Feb 1991, pp 131--133.
     */
    slong n;
    int invertible;
    fmpq_mat_t S, I, A, B;
    if (!fmpq_mat_is_square(Q)) abort(); /* assert */
    n = fmpq_mat_nrows(Q);
    fmpq_mat_init(S, n, n);
    fmpq_mat_init(I, n, n);
    fmpq_mat_init(A, n, n);
    fmpq_mat_init(B, n, n);
    _fmpq_mat_randtest_skew_symmetric(S, state, bits);
    fmpq_mat_one(I);
    fmpq_mat_sub(A, S, I);
    fmpq_mat_add(B, S, I);
    invertible = fmpq_mat_solve_fraction_free(Q, A, B);
    if (!invertible) abort(); /* assert */
    fmpq_mat_clear(S);
    fmpq_mat_clear(I);
    fmpq_mat_clear(A);
    fmpq_mat_clear(B);
}

void
_test_eigenvalue_separation(flint_rand_t state)
{
    slong i, n, prec, k;
    fmpq_mat_t Q, QT, D, A, Dq;
    mp_bitcnt_t bits;
    int ambiguous;

    bits = 20;
    n = 4;

    fmpq_mat_init(Q, n, n);
    fmpq_mat_init(QT, n, n);
    fmpq_mat_init(A, n, n);
    fmpq_mat_init(D, n, n);
    fmpq_mat_init(Dq, n, n);

    _fmpq_mat_randtest_orthogonal(Q, state, bits);
    fmpq_mat_transpose(QT, Q);

    /* check an orthogonality property */
    fmpq_mat_mul(A, Q, QT);
    if (!_fmpq_mat_is_one(A)) abort(); /* assert */

    k = 1 << 24;
    fmpq_set_si(fmpq_mat_entry(D, 0, 0), 1, k);
    fmpq_set_si(fmpq_mat_entry(D, 1, 1), k, k);
    fmpq_set_si(fmpq_mat_entry(D, 2, 2), k+1, k);
    fmpq_set_si(fmpq_mat_entry(D, 3, 3), k, 1);
    fmpq_mat_mul(A, Q, D);
    fmpq_mat_mul(A, A, QT);

    if (!_fmpq_mat_is_symmetric(A)) abort(); /* assert */

    fmpq_mat_init(Dq, n, 1);
    for (i = 0; i < n; i++)
    {
        fmpq_set(fmpq_mat_entry(Dq, i, 0), fmpq_mat_entry(D, i, i));
    }

    /* try to separate the eigenvalues in A */
    ambiguous = 0;
    for (prec = 2; prec < 10000; prec <<= 1)
    {
        arb_mat_t Ar, Dr, Pr;

        /* flint_printf("prec=%wd\n", prec); */
        arb_mat_init(Ar, n, n);
        arb_mat_init(Dr, n, 1);
        arb_mat_init(Pr, n, n);
        arb_mat_set_fmpq_mat(Ar, A, prec);
        /* arb_mat_printd(Ar, 15); */
        ambiguous = arb_mat_symmetric_diagonalization(Dr, Pr, Ar, prec);
        if (!ambiguous)
        {
            if (!arb_mat_contains_fmpq_mat(Dr, Dq))
            {
                flint_printf("FAIL (eigenvalue containment)\n");
                flint_printf("Dr = \n"); arb_mat_printd(Dr, 15); flint_printf("\n");
                flint_printf("Dq = \n"); fmpq_mat_print(Dq); flint_printf("\n");
                abort();
            }
        }
        arb_mat_clear(Ar);
        arb_mat_clear(Dr);
        arb_mat_clear(Pr);
        if (!ambiguous)
            break;
    }
    if (ambiguous)
    {
        flint_printf("FAIL (eigenvalue ambiguity)\n");
        flint_printf("Q = \n"); fmpq_mat_print(Q); flint_printf("\n");
        abort();
    }

    fmpq_mat_clear(Q);
    fmpq_mat_clear(QT);
    fmpq_mat_clear(A);
    fmpq_mat_clear(D);
    fmpq_mat_clear(Dq);
}

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("symmetric_diagonalization....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100; iter++)
    {
        _test_eigenvalue_separation(state);
    }

    /* test known diagonalizations */
    for (iter = 0; iter < 10000; iter++)
    {
        _test_dct2(state);
    }

    for (iter = 0; iter < 10000; iter++)
    {
        slong i, j, n;
        slong prec;
        arb_mat_t D, P, A;

        n = n_randint(state, 8);
        prec = 2 + n_randint(state, 202);

        arb_mat_init(D, n, 1);
        arb_mat_init(P, n, n);
        arb_mat_init(A, n, n);

        arb_mat_randtest(A, state, 2 + n_randint(state, 100), 10);
        for (i = 0; i < n; i++)
            for (j = 0; j < i; j++)
                arb_set(arb_mat_entry(A, i, j), arb_mat_entry(A, j, i));

        arb_mat_symmetric_diagonalization(D, P, A, prec);

        /* charpoly should contain zero at eigenvalues of A */
        {
            slong i;
            arb_poly_t f;
            arb_t y;

            arb_poly_init(f);
            arb_init(y);

            arb_mat_charpoly(f, A, prec);

            for (i = 0; i < n; i++)
            {
                arb_srcptr x = arb_mat_entry(D, i, 0);

                arb_poly_evaluate(y, f, x, prec);
                if (!arb_contains_zero(y))
                {
                    flint_printf("FAIL (charpoly)\n");
                    flint_printf("A = \n"); arb_mat_printd(A, 15); flint_printf("\n");
                    flint_printf("D = \n"); arb_mat_printd(D, 15); flint_printf("\n");
                    flint_printf("P = \n"); arb_mat_printd(P, 15); flint_printf("\n");
                    abort();
                }
            }

            arb_poly_clear(f);
            arb_clear(y);
        }

        /* Multiplying out the decomposition should give the original matrix */
        {
            /* A = P * D * P^T */
            slong i, j, k;
            arb_t t;
            arb_mat_t Y;

            arb_init(t);
            arb_mat_init(Y, n, n);

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    for (k = 0; k < n; k++)
                    {
                        arb_mul(t, arb_mat_entry(P, i, k),
                                   arb_mat_entry(P, j, k), prec);
                        arb_addmul(arb_mat_entry(Y, i, j),
                                   arb_mat_entry(D, k, 0), t, prec);
                    }
                }
            }

            if (!arb_mat_contains(Y, A))
            {
                flint_printf("FAIL (decomposition expansion)\n");
                flint_printf("A = \n"); arb_mat_printd(A, 15);
                flint_printf("D = \n"); arb_mat_printd(D, 15);
                flint_printf("P = \n"); arb_mat_printd(P, 15);
                flint_printf("P * D * P^T = \n"); arb_mat_printd(Y, 15);
                abort();
            }

            arb_clear(t);
            arb_mat_clear(Y);
        }

        arb_mat_clear(D);
        arb_mat_clear(P);
        arb_mat_clear(A);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
