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

#ifndef ARB_MAT_H
#define ARB_MAT_H

#ifdef ARB_MAT_INLINES_C
#define ARB_MAT_INLINE
#else
#define ARB_MAT_INLINE static __inline__
#endif

#include "arb.h"
#include "fmpz_mat.h"
#include "fmpq_mat.h"
#include "perm.h"
#include "arb_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    arb_ptr entries;
    long r;
    long c;
    arb_ptr * rows;
}
arb_mat_struct;

typedef arb_mat_struct arb_mat_t[1];

#define arb_mat_entry(mat,i,j) ((mat)->rows[i] + (j))
#define arb_mat_nrows(mat) ((mat)->r)
#define arb_mat_ncols(mat) ((mat)->c)


/* Memory management */

void arb_mat_init(arb_mat_t mat, long r, long c);

void arb_mat_clear(arb_mat_t mat);

ARB_MAT_INLINE void
arb_mat_swap(arb_mat_t mat1, arb_mat_t mat2)
{
    arb_mat_struct t = *mat1;
    *mat1 = *mat2;
    *mat2 = t;
}

/* Conversions */

void arb_mat_set(arb_mat_t dest, const arb_mat_t src);

void arb_mat_set_fmpz_mat(arb_mat_t dest, const fmpz_mat_t src);

void arb_mat_set_fmpq_mat(arb_mat_t dest, const fmpq_mat_t src, long prec);

/* Random generation */

void arb_mat_randtest(arb_mat_t mat, flint_rand_t state, long prec, long mag_bits);

/* I/O */

void arb_mat_printd(const arb_mat_t mat, long digits);

/* Comparisons */

int arb_mat_equal(const arb_mat_t mat1, const arb_mat_t mat2);

int arb_mat_overlaps(const arb_mat_t mat1, const arb_mat_t mat2);

int arb_mat_contains(const arb_mat_t mat1, const arb_mat_t mat2);

int arb_mat_contains_fmpq_mat(const arb_mat_t mat1, const fmpq_mat_t mat2);

int arb_mat_contains_fmpz_mat(const arb_mat_t mat1, const fmpz_mat_t mat2);

/* Special matrices */

void arb_mat_zero(arb_mat_t mat);

void arb_mat_one(arb_mat_t mat);

/* Norms */

void arb_mat_bound_inf_norm(mag_t b, const arb_mat_t A);

/* Arithmetic */

void arb_mat_neg(arb_mat_t dest, const arb_mat_t src);

void arb_mat_add(arb_mat_t res, const arb_mat_t mat1, const arb_mat_t mat2, long prec);

void arb_mat_sub(arb_mat_t res, const arb_mat_t mat1, const arb_mat_t mat2, long prec);

void arb_mat_mul(arb_mat_t res, const arb_mat_t mat1, const arb_mat_t mat2, long prec);

void arb_mat_mul_classical(arb_mat_t C, const arb_mat_t A, const arb_mat_t B, long prec);

void arb_mat_mul_threaded(arb_mat_t C, const arb_mat_t A, const arb_mat_t B, long prec);

void arb_mat_pow_ui(arb_mat_t B, const arb_mat_t A, ulong exp, long prec);

/* Scalar arithmetic */

ARB_MAT_INLINE void
arb_mat_scalar_mul_2exp_si(arb_mat_t B, const arb_mat_t A, long c)
{
    long i, j;

    for (i = 0; i < arb_mat_nrows(A); i++)
        for (j = 0; j < arb_mat_ncols(A); j++)
            arb_mul_2exp_si(arb_mat_entry(B, i, j), arb_mat_entry(A, i, j), c);
}

ARB_MAT_INLINE void
arb_mat_scalar_addmul_si(arb_mat_t B, const arb_mat_t A, long c, long prec)
{
    long i, j;

    for (i = 0; i < arb_mat_nrows(A); i++)
        for (j = 0; j < arb_mat_ncols(A); j++)
            arb_addmul_si(arb_mat_entry(B, i, j), arb_mat_entry(A, i, j), c, prec);
}

ARB_MAT_INLINE void
arb_mat_scalar_mul_si(arb_mat_t B, const arb_mat_t A, long c, long prec)
{
    long i, j;

    for (i = 0; i < arb_mat_nrows(A); i++)
        for (j = 0; j < arb_mat_ncols(A); j++)
            arb_mul_si(arb_mat_entry(B, i, j), arb_mat_entry(A, i, j), c, prec);
}

ARB_MAT_INLINE void
arb_mat_scalar_div_si(arb_mat_t B, const arb_mat_t A, long c, long prec)
{
    long i, j;

    for (i = 0; i < arb_mat_nrows(A); i++)
        for (j = 0; j < arb_mat_ncols(A); j++)
            arb_div_si(arb_mat_entry(B, i, j), arb_mat_entry(A, i, j), c, prec);
}

ARB_MAT_INLINE void
arb_mat_scalar_addmul_fmpz(arb_mat_t B, const arb_mat_t A, const fmpz_t c, long prec)
{
    long i, j;

    for (i = 0; i < arb_mat_nrows(A); i++)
        for (j = 0; j < arb_mat_ncols(A); j++)
            arb_addmul_fmpz(arb_mat_entry(B, i, j), arb_mat_entry(A, i, j), c, prec);
}

ARB_MAT_INLINE void
arb_mat_scalar_mul_fmpz(arb_mat_t B, const arb_mat_t A, const fmpz_t c, long prec)
{
    long i, j;

    for (i = 0; i < arb_mat_nrows(A); i++)
        for (j = 0; j < arb_mat_ncols(A); j++)
            arb_mul_fmpz(arb_mat_entry(B, i, j), arb_mat_entry(A, i, j), c, prec);
}

ARB_MAT_INLINE void
arb_mat_scalar_div_fmpz(arb_mat_t B, const arb_mat_t A, const fmpz_t c, long prec)
{
    long i, j;

    for (i = 0; i < arb_mat_nrows(A); i++)
        for (j = 0; j < arb_mat_ncols(A); j++)
            arb_div_fmpz(arb_mat_entry(B, i, j), arb_mat_entry(A, i, j), c, prec);
}

ARB_MAT_INLINE void
arb_mat_scalar_addmul_arb(arb_mat_t B, const arb_mat_t A, const arb_t c, long prec)
{
    long i, j;

    for (i = 0; i < arb_mat_nrows(A); i++)
        for (j = 0; j < arb_mat_ncols(A); j++)
            arb_addmul(arb_mat_entry(B, i, j), arb_mat_entry(A, i, j), c, prec);
}

ARB_MAT_INLINE void
arb_mat_scalar_mul_arb(arb_mat_t B, const arb_mat_t A, const arb_t c, long prec)
{
    long i, j;

    for (i = 0; i < arb_mat_nrows(A); i++)
        for (j = 0; j < arb_mat_ncols(A); j++)
            arb_mul(arb_mat_entry(B, i, j), arb_mat_entry(A, i, j), c, prec);
}

ARB_MAT_INLINE void
arb_mat_scalar_div_arb(arb_mat_t B, const arb_mat_t A, const arb_t c, long prec)
{
    long i, j;

    for (i = 0; i < arb_mat_nrows(A); i++)
        for (j = 0; j < arb_mat_ncols(A); j++)
            arb_div(arb_mat_entry(B, i, j), arb_mat_entry(A, i, j), c, prec);
}

/* Solving */

ARB_MAT_INLINE void
arb_mat_swap_rows(arb_mat_t mat, long * perm, long r, long s)
{
    if (r != s)
    {
        arb_ptr u;
        long t;

        if (perm != NULL)
        {
            t = perm[s];
            perm[s] = perm[r];
            perm[r] = t;
        }

        u = mat->rows[s];
        mat->rows[s] = mat->rows[r];
        mat->rows[r] = u;
    }
}

long arb_mat_find_pivot_partial(const arb_mat_t mat,
                                    long start_row, long end_row, long c);

int arb_mat_lu(long * P, arb_mat_t LU, const arb_mat_t A, long prec);

void arb_mat_solve_lu_precomp(arb_mat_t X, const long * perm,
    const arb_mat_t A, const arb_mat_t B, long prec);

int arb_mat_solve(arb_mat_t X, const arb_mat_t A, const arb_mat_t B, long prec);

int arb_mat_inv(arb_mat_t X, const arb_mat_t A, long prec);

void arb_mat_det(arb_t det, const arb_mat_t A, long prec);

/* Special functions */

void arb_mat_exp(arb_mat_t B, const arb_mat_t A, long prec);

void _arb_mat_charpoly(arb_ptr cp, const arb_mat_t mat, long prec);

void arb_mat_charpoly(arb_poly_t cp, const arb_mat_t mat, long prec);

#ifdef __cplusplus
}
#endif

#endif

