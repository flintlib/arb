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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#ifndef PARTITIONS_H
#define PARTITIONS_H

#include <math.h>
#include "flint.h"
#include "arith.h"
#include "arb.h"

#ifdef __cplusplus
extern "C" {
#endif

void partitions_rademacher_bound(arf_t b, const fmpz_t n, ulong N);

void partitions_hrr_sum_arb(arb_t x, const fmpz_t n, long N0, long N, int use_doubles);

void partitions_fmpz_fmpz(fmpz_t p, const fmpz_t n, int use_doubles);

void partitions_fmpz_ui(fmpz_t p, ulong n);

void partitions_fmpz_ui_using_doubles(fmpz_t p, ulong n);

#ifdef __cplusplus
}
#endif

#endif
