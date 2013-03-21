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

#ifndef ZETA_H
#define ZETA_H

#include <math.h>
#include "flint.h"
#include "fmprb.h"
#include "fmpcb.h"

void fmprb_const_zeta3_bsplit(fmprb_t x, long prec);

void fmprb_const_khinchin(fmprb_t K, long prec);

void fmprb_zeta_ui_asymp(fmprb_t x, ulong s, long prec);
void fmprb_zeta_ui_bsplit(fmprb_t x, ulong s, long prec);
void fmprb_zeta_ui_euler_product(fmprb_t z, ulong s, long prec);
void fmprb_zeta_ui_bernoulli(fmprb_t x, ulong n, long prec);
void fmprb_zeta_ui_vec_borwein(fmprb_struct * z, ulong start, long num, ulong step, long prec);
void fmprb_zeta_ui(fmprb_t x, ulong n, long prec);

void fmprb_zeta_ui_vec_even(fmprb_struct * x, ulong start, long num, long prec);
void fmprb_zeta_ui_vec_odd(fmprb_struct * x, ulong start, long num, long prec);
void fmprb_zeta_ui_vec(fmprb_struct * x, ulong start, long num, long prec);

void fmpcb_zeta_series_em_sum(fmpcb_struct * z, const fmpcb_t s, const fmpcb_t a, int deflate, ulong N, ulong M, long d, long prec);
void fmpcb_zeta_series_em_choose_param(fmpr_t bound, ulong * N, ulong * M, const fmpcb_t s, const fmpcb_t a, long d, long target, long prec);
void fmpcb_zeta_series_em_bound(fmpr_t bound, const fmpcb_t s, const fmpcb_t a, long N, long M, long d, long wp);
void fmpcb_zeta_series_em_vec_bound(fmprb_struct * vec, const fmpcb_t s, const fmpcb_t a, ulong N, ulong M, long d, long wp);
void fmpcb_zeta_series(fmpcb_struct * z, const fmpcb_t s, const fmpcb_t a, int deflate, long d, long prec);
void fmpcb_zeta(fmpcb_t z, const fmpcb_t s, long prec);

#endif

