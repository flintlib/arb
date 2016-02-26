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

#ifndef FMPZ_MAT_EXTRAS_H
#define FMPZ_MAT_EXTRAS_H

#include "flint.h"
#include "fmpz_mat.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Convenience functions related to sparsity structure */

int fmpz_mat_is_hollow(const fmpz_mat_t mat);

int fmpz_mat_is_diagonal(const fmpz_mat_t mat);

int fmpz_mat_is_nonnegative(const fmpz_mat_t mat);

int fmpz_mat_is_lower_triangular(const fmpz_mat_t mat);

void fmpz_mat_entrywise_not_is_zero(fmpz_mat_t dest, const fmpz_mat_t src);

slong fmpz_mat_count_nonzero(const fmpz_mat_t mat);

/* Arithmetic */

void fmpz_mat_add_ui_entrywise(fmpz_mat_t B, const fmpz_mat_t A, ulong x);

void fmpz_mat_sub_ui_entrywise(fmpz_mat_t B, const fmpz_mat_t A, ulong x);

/* Graph theory */

void fmpz_mat_transitive_closure(fmpz_mat_t B, const fmpz_mat_t A);

void fmpz_mat_unweighted_all_pairs_longest_walk(fmpz_mat_t B, const fmpz_mat_t A);

#ifdef __cplusplus
}
#endif

#endif

