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

void fmpz_mat_transitive_closure(fmpz_mat_t dest, const fmpz_mat_t src);

void fmpz_mat_entrywise_nilpotence_degree(fmpz_mat_t dest, const fmpz_mat_t src);

void fmpz_mat_entrywise_not_is_zero(fmpz_mat_t dest, const fmpz_mat_t src);

int fmpz_mat_is_hollow(const fmpz_mat_t mat);

int fmpz_mat_is_nonnegative(const fmpz_mat_t mat);

int fmpz_mat_is_lower_triangular(const fmpz_mat_t mat);

#ifdef __cplusplus
}
#endif

#endif

