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

#ifndef FMPCB_MAT_H
#define FMPCB_MAT_H

#include "fmprb.h"
#include "fmpcb.h"
#include "fmpz_mat.h"
#include "fmpq_mat.h"
#include "fmprb_mat.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct
{
    fmpcb_ptr entries;
    long r;
    long c;
    fmpcb_ptr * rows;
}
fmpcb_mat_struct;

typedef fmpcb_mat_struct fmpcb_mat_t[1];

#define fmpcb_mat_entry(mat,i,j) ((mat)->rows[i] + (j))
#define fmpcb_mat_nrows(mat) ((mat)->r)
#define fmpcb_mat_ncols(mat) ((mat)->c)


/* Memory management */

void fmpcb_mat_init(fmpcb_mat_t mat, long r, long c);

void fmpcb_mat_clear(fmpcb_mat_t mat);

static __inline__ void
fmpcb_mat_swap(fmpcb_mat_t mat1, fmpcb_mat_t mat2)
{
    fmpcb_mat_struct t = *mat1;
    *mat1 = *mat2;
    *mat2 = t;
}

/* Conversions */

void fmpcb_mat_set(fmpcb_mat_t dest, const fmpcb_mat_t src);

void fmpcb_mat_set_fmpz_mat(fmpcb_mat_t dest, const fmpz_mat_t src);

void fmpcb_mat_set_fmpq_mat(fmpcb_mat_t dest, const fmpq_mat_t src, long prec);

/* I/O */

void fmpcb_mat_printd(const fmpcb_mat_t mat, long digits);

/* Comparisons */

int fmpcb_mat_equal(const fmpcb_mat_t mat1, const fmpcb_mat_t mat2);

int fmpcb_mat_overlaps(const fmpcb_mat_t mat1, const fmpcb_mat_t mat2);

int fmpcb_mat_contains(const fmpcb_mat_t mat1, const fmpcb_mat_t mat2);

int fmpcb_mat_contains_fmpq_mat(const fmpcb_mat_t mat1, const fmpq_mat_t mat2);

int fmpcb_mat_contains_fmpz_mat(const fmpcb_mat_t mat1, const fmpz_mat_t mat2);

/* Special matrices */

void fmpcb_mat_zero(fmpcb_mat_t mat);

void fmpcb_mat_one(fmpcb_mat_t mat);

/* Arithmetic */

void fmpcb_mat_neg(fmpcb_mat_t dest, const fmpcb_mat_t src);

void fmpcb_mat_add(fmpcb_mat_t res, const fmpcb_mat_t mat1, const fmpcb_mat_t mat2, long prec);

void fmpcb_mat_sub(fmpcb_mat_t res, const fmpcb_mat_t mat1, const fmpcb_mat_t mat2, long prec);

void fmpcb_mat_mul(fmpcb_mat_t res, const fmpcb_mat_t mat1, const fmpcb_mat_t mat2, long prec);

void fmpcb_mat_pow_ui(fmpcb_mat_t B, const fmpcb_mat_t A, ulong exp, long prec);

/* Solving */

static __inline__ void
fmpcb_mat_swap_rows(fmpcb_mat_t mat, long * perm, long r, long s)
{
    if (r != s)
    {
        fmpcb_ptr u;
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

long fmpcb_mat_find_pivot_partial(const fmpcb_mat_t mat,
                                    long start_row, long end_row, long c);

int fmpcb_mat_lu(long * P, fmpcb_mat_t LU, const fmpcb_mat_t A, long prec);

void fmpcb_mat_solve_lu_precomp(fmpcb_mat_t X, const long * perm,
    const fmpcb_mat_t A, const fmpcb_mat_t B, long prec);

int fmpcb_mat_solve(fmpcb_mat_t X, const fmpcb_mat_t A, const fmpcb_mat_t B, long prec);

int fmpcb_mat_inv(fmpcb_mat_t X, const fmpcb_mat_t A, long prec);

void fmpcb_mat_det(fmpcb_t det, const fmpcb_mat_t A, long prec);


#ifdef __cplusplus
}
#endif

#endif

