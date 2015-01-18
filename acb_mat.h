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

#ifndef ACB_MAT_H
#define ACB_MAT_H

#ifdef ACB_MAT_INLINES_C
#define ACB_MAT_INLINE
#else
#define ACB_MAT_INLINE static __inline__
#endif

#include "arb.h"
#include "acb.h"
#include "fmpz_mat.h"
#include "fmpq_mat.h"
#include "arb_mat.h"
#include "acb_poly.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct
{
    acb_ptr entries;
    long r;
    long c;
    acb_ptr * rows;
}
acb_mat_struct;

typedef acb_mat_struct acb_mat_t[1];

#define acb_mat_entry(mat,i,j) ((mat)->rows[i] + (j))
#define acb_mat_nrows(mat) ((mat)->r)
#define acb_mat_ncols(mat) ((mat)->c)


/* Memory management */

void acb_mat_init(acb_mat_t mat, long r, long c);

void acb_mat_clear(acb_mat_t mat);

ACB_MAT_INLINE void
acb_mat_swap(acb_mat_t mat1, acb_mat_t mat2)
{
    acb_mat_struct t = *mat1;
    *mat1 = *mat2;
    *mat2 = t;
}

/* Conversions */

void acb_mat_set(acb_mat_t dest, const acb_mat_t src);

void acb_mat_set_fmpz_mat(acb_mat_t dest, const fmpz_mat_t src);

void acb_mat_set_fmpq_mat(acb_mat_t dest, const fmpq_mat_t src, long prec);

/* Random generation */

void acb_mat_randtest(acb_mat_t mat, flint_rand_t state, long prec, long mag_bits);

/* I/O */

void acb_mat_printd(const acb_mat_t mat, long digits);

/* Comparisons */

int acb_mat_equal(const acb_mat_t mat1, const acb_mat_t mat2);

int acb_mat_overlaps(const acb_mat_t mat1, const acb_mat_t mat2);

int acb_mat_contains(const acb_mat_t mat1, const acb_mat_t mat2);

int acb_mat_contains_fmpq_mat(const acb_mat_t mat1, const fmpq_mat_t mat2);

int acb_mat_contains_fmpz_mat(const acb_mat_t mat1, const fmpz_mat_t mat2);

int acb_mat_is_real(const acb_mat_t mat);

/* Special matrices */

void acb_mat_zero(acb_mat_t mat);

void acb_mat_one(acb_mat_t mat);

/* Norms */

void acb_mat_bound_inf_norm(mag_t b, const acb_mat_t A);

/* Arithmetic */

void acb_mat_neg(acb_mat_t dest, const acb_mat_t src);

void acb_mat_add(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, long prec);

void acb_mat_sub(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, long prec);

void acb_mat_mul(acb_mat_t res, const acb_mat_t mat1, const acb_mat_t mat2, long prec);

void acb_mat_pow_ui(acb_mat_t B, const acb_mat_t A, ulong exp, long prec);

/* Scalar arithmetic */

ACB_MAT_INLINE void
acb_mat_scalar_mul_2exp_si(acb_mat_t B, const acb_mat_t A, long c)
{
    long i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            acb_mul_2exp_si(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), c);
}

ACB_MAT_INLINE void
acb_mat_scalar_addmul_si(acb_mat_t B, const acb_mat_t A, long c, long prec)
{
    long i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            acb_addmul_si(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), c, prec);
}

ACB_MAT_INLINE void
acb_mat_scalar_mul_si(acb_mat_t B, const acb_mat_t A, long c, long prec)
{
    long i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            acb_mul_si(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), c, prec);
}

ACB_MAT_INLINE void
acb_mat_scalar_div_si(acb_mat_t B, const acb_mat_t A, long c, long prec)
{
    long i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            acb_div_si(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), c, prec);
}

ACB_MAT_INLINE void
acb_mat_scalar_addmul_fmpz(acb_mat_t B, const acb_mat_t A, const fmpz_t c, long prec)
{
    long i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            acb_addmul_fmpz(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), c, prec);
}

ACB_MAT_INLINE void
acb_mat_scalar_mul_fmpz(acb_mat_t B, const acb_mat_t A, const fmpz_t c, long prec)
{
    long i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            acb_mul_fmpz(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), c, prec);
}

ACB_MAT_INLINE void
acb_mat_scalar_div_fmpz(acb_mat_t B, const acb_mat_t A, const fmpz_t c, long prec)
{
    long i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            acb_div_fmpz(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), c, prec);
}

ACB_MAT_INLINE void
acb_mat_scalar_addmul_acb(acb_mat_t B, const acb_mat_t A, const acb_t c, long prec)
{
    long i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            acb_addmul(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), c, prec);
}

ACB_MAT_INLINE void
acb_mat_scalar_mul_acb(acb_mat_t B, const acb_mat_t A, const acb_t c, long prec)
{
    long i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            acb_mul(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), c, prec);
}

ACB_MAT_INLINE void
acb_mat_scalar_div_acb(acb_mat_t B, const acb_mat_t A, const acb_t c, long prec)
{
    long i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            acb_div(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), c, prec);
}

ACB_MAT_INLINE void
acb_mat_scalar_addmul_arb(acb_mat_t B, const acb_mat_t A, const arb_t c, long prec)
{
    long i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            acb_addmul_arb(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), c, prec);
}

ACB_MAT_INLINE void
acb_mat_scalar_mul_arb(acb_mat_t B, const acb_mat_t A, const arb_t c, long prec)
{
    long i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            acb_mul_arb(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), c, prec);
}

ACB_MAT_INLINE void
acb_mat_scalar_div_arb(acb_mat_t B, const acb_mat_t A, const arb_t c, long prec)
{
    long i, j;

    for (i = 0; i < acb_mat_nrows(A); i++)
        for (j = 0; j < acb_mat_ncols(A); j++)
            acb_div_arb(acb_mat_entry(B, i, j), acb_mat_entry(A, i, j), c, prec);
}

/* Solving */

ACB_MAT_INLINE void
acb_mat_swap_rows(acb_mat_t mat, long * perm, long r, long s)
{
    if (r != s)
    {
        acb_ptr u;
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

long acb_mat_find_pivot_partial(const acb_mat_t mat,
                                    long start_row, long end_row, long c);

int acb_mat_lu(long * P, acb_mat_t LU, const acb_mat_t A, long prec);

void acb_mat_solve_lu_precomp(acb_mat_t X, const long * perm,
    const acb_mat_t A, const acb_mat_t B, long prec);

int acb_mat_solve(acb_mat_t X, const acb_mat_t A, const acb_mat_t B, long prec);

int acb_mat_inv(acb_mat_t X, const acb_mat_t A, long prec);

void acb_mat_det(acb_t det, const acb_mat_t A, long prec);

/* Special functions */

void acb_mat_exp(acb_mat_t B, const acb_mat_t A, long prec);

void _acb_mat_charpoly(acb_ptr cp, const acb_mat_t mat, long prec);

void acb_mat_charpoly(acb_poly_t cp, const acb_mat_t mat, long prec);

#ifdef __cplusplus
}
#endif

#endif

