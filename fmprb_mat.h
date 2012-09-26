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

#ifndef FMPRB_MAT_H
#define FMPRB_MAT_H

#include "fmprb.h"
#include "fmpz_mat.h"
#include "fmpq_mat.h"

typedef struct
{
    fmprb_struct * entries;
    long r;
    long c;
    fmprb_struct ** rows;
}
fmprb_mat_struct;

typedef fmprb_mat_struct fmprb_mat_t[1];

#define fmprb_mat_entry(mat,i,j) ((mat)->rows[i] + (j))
#define fmprb_mat_nrows(mat) ((mat)->r)
#define fmprb_mat_ncols(mat) ((mat)->c)


/* Memory management */

void fmprb_mat_init(fmprb_mat_t mat, long r, long c);

void fmprb_mat_clear(fmprb_mat_t mat);

/* Conversions */

void fmprb_mat_set(fmprb_mat_t dest, const fmprb_mat_t src);

void fmprb_mat_set_fmpz_mat(fmprb_mat_t dest, const fmpz_mat_t src);

void fmprb_mat_set_fmpq_mat(fmprb_mat_t dest, const fmpq_mat_t src, long prec);

/* I/O */

void fmprb_mat_printd(const fmprb_mat_t mat, long digits);

/* Comparisons */

int fmprb_mat_equal(const fmprb_mat_t mat1, const fmprb_mat_t mat2);

int fmprb_mat_contains_fmpq_mat(const fmprb_mat_t mat1, const fmpq_mat_t mat2);

/* Special matrices */

void fmprb_mat_zero(fmprb_mat_t mat);

void fmprb_mat_one(fmprb_mat_t mat);

/* Arithmetic */

void fmprb_mat_neg(fmprb_mat_t dest, const fmprb_mat_t src);

void fmprb_mat_add(fmprb_mat_t res, const fmprb_mat_t mat1, const fmprb_mat_t mat2, long prec);

void fmprb_mat_sub(fmprb_mat_t res, const fmprb_mat_t mat1, const fmprb_mat_t mat2, long prec);

void fmprb_mat_mul(fmprb_mat_t res, const fmprb_mat_t mat1, const fmprb_mat_t mat2, long prec);

#endif

