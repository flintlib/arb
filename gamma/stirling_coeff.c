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

#include "gamma.h"
#include "bernoulli.h"

void
gamma_stirling_coeff(fmprb_t b, ulong k, int digamma, long prec)
{
    fmpz_t d;
    fmpz_init(d);
    BERNOULLI_ENSURE_CACHED(2 * k);
    fmprb_set_round_fmpz(b, fmpq_numref(bernoulli_cache + 2 * k), prec);
    if (digamma)
        fmpz_mul_ui(d, fmpq_denref(bernoulli_cache + 2 * k), 2 * k);
    else
        fmpz_mul2_uiui(d, fmpq_denref(bernoulli_cache + 2 * k), 2 * k, 2 * k - 1);
    fmprb_div_fmpz(b, b, d, prec);
    fmpz_clear(d);
}

