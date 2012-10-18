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

#include "fmpz_holonomic.h"
#include "nmod_poly_mat.h"

static void
mat_bsplit(nmod_poly_mat_t M, nmod_poly_t Q, nmod_poly_mat_t C, nmod_poly_t Cden, ulong a, ulong b)
{
    long i, j, r = C->r;

    if (b - a == 1)
    {
        for (i = 0; i < r; i++)
            for (j = 0; j < r; j++)
                nmod_poly_taylor_shift(M->rows[i] + j, C->rows[i] + j, a);

        nmod_poly_taylor_shift(Q, Cden, a);
    }
    else
    {
        nmod_poly_mat_t M1, M2;
        nmod_poly_t Q1, Q2;
        ulong m = a + (b - a) / 2;

        nmod_poly_mat_init(M1, r, r, Q->mod.n);
        nmod_poly_mat_init(M2, r, r, Q->mod.n);

        nmod_poly_init(Q1, Q->mod.n);
        nmod_poly_init(Q2, Q->mod.n);

        mat_bsplit(M1, Q1, C, Cden, a, m);
        mat_bsplit(M2, Q2, C, Cden, m, b);

        nmod_poly_mat_mul(M, M2, M1);
        nmod_poly_mul(Q, Q2, Q1);

        nmod_poly_clear(Q1);
        nmod_poly_clear(Q2);

        nmod_poly_mat_clear(M1);
        nmod_poly_mat_clear(M2);
    }
}

static __inline__ void
nmod_mat_swap(nmod_mat_t A, nmod_mat_t B)
{
    nmod_mat_struct T = *A;
    *A = *B;
    *B = T;
}

static __inline__ void
nmod_mat_one(nmod_mat_t A)
{
    long i, j;

    for (i = 0; i < A->r; i++)
        for (j = 0; j < A->c; j++)
            A->rows[i][j] = (i == j);
}


void
fmpz_holonomic_forward_nmod_mat(nmod_mat_t M, mp_limb_t * Q, const fmpz_holonomic_t op, ulong start, ulong n)
{
    long i, j, k, r;
    nmod_poly_mat_t C;
    nmod_poly_mat_t PM;
    nmod_poly_t PQ, Cden;
    nmod_mat_t T, U, *Y, *X;
    mp_ptr ys, xs;
    mp_limb_t p, q;
    ulong m;

    nmod_mat_one(M);
    q = 1;

    r = fmpz_holonomic_order(op);
    p = M->mod.n;

    /* number of evaluation points */
    m = n_sqrt(n);

    nmod_poly_mat_init(C, r, r, p);
    nmod_poly_mat_init(PM, r, r, p);
    nmod_poly_init(PQ, p);
    nmod_poly_init(Cden, p);
    nmod_mat_init(T, r, r, p);
    nmod_mat_init(U, r, r, p);

    /* construct companion matrix and denominator */
    fmpz_poly_get_nmod_poly(Cden, op->coeffs + r);
    nmod_poly_neg(Cden, Cden);
    for (i = 0; i < r - 1; i++)
        nmod_poly_set(C->rows[i] + i + 1, Cden);
    for (i = 0; i < r; i++)
        fmpz_poly_get_nmod_poly(C->rows[r - 1] + i, op->coeffs + i);

    if (start != 0)
    {
        for (i = 0; i < r; i++)
            for (j = 0; j < r; j++)
                nmod_poly_taylor_shift(C->rows[i] + j, C->rows[i] + j, start);
        nmod_poly_taylor_shift(Cden, Cden, start);
    }

    if (m > 0)
    {
        xs = _nmod_vec_init(m);
        ys = _nmod_vec_init(m);

        X = flint_malloc(sizeof(nmod_mat_t) * m);
        Y = flint_malloc(sizeof(nmod_mat_t) * m);

        for (k = 0; k < m; k++)
        {
            nmod_mat_init(X[k], r, r, p);
            nmod_mat_init(Y[k], r, r, p);
        }

        /* compute product of companion matrices */
        mat_bsplit(PM, PQ, C, Cden, 0, m);

        /* points for multipoint evaluation */
        for (k = 0; k < m; k++)
            xs[k] = (k * m) % p;

        /* multipoint evaluation of numerator */
        for (i = 0; i < r; i++)
        {
            for (j = 0; j < r; j++)
            {
                nmod_poly_evaluate_nmod_vec(ys, PM->rows[i] + j, xs, m);
                for (k = 0; k < m; k++)
                    X[k]->rows[i][j] = ys[k];
            }
        }

        /* multipoint evaluation of denominator */
        nmod_poly_evaluate_nmod_vec(ys, PQ, xs, m);
        q = 1;
        for (k = 0; k < m; k++)
            q = n_mulmod2_preinv(q, ys[k], M->mod.n, M->mod.ninv);

        /* multiply together evaluated matrices */
        nmod_mat_set(M, X[0]);
        for (i = 1; i < m; i++)
        {
            nmod_mat_mul(T, X[i], M);
            nmod_mat_swap(M, T);
        }

        _nmod_vec_clear(xs);
        _nmod_vec_clear(ys);

        for (k = 0; k < m; k++)
        {
            nmod_mat_clear(X[k]);
            nmod_mat_clear(Y[k]);
        }

        flint_free(X);
        flint_free(Y);
    }

    /* fill in the rest */
    for (i = m * m; i < n; i++)
    {
        nmod_poly_mat_evaluate_nmod(U, C, i % p);
        nmod_mat_mul(T, U, M);
        nmod_mat_swap(M, T);
        q = n_mulmod2_preinv(q,
            nmod_poly_evaluate_nmod(Cden, i % p), M->mod.n, M->mod.ninv);
    }

    *Q = q;

    nmod_poly_mat_clear(C);
    nmod_poly_mat_clear(PM);
    nmod_poly_clear(PQ);
    nmod_poly_clear(Cden);
    nmod_mat_clear(T);
    nmod_mat_clear(U);
}

