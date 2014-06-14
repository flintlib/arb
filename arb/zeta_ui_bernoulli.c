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

#include "arb.h"
#include "bernoulli.h"

void
arb_zeta_ui_bernoulli(arb_t x, ulong n, long prec)
{
    fmpq_t b;
    arb_t t, f;
    long wp;

    if (n % 2)
        abort();

    wp = prec + FLINT_BIT_COUNT(n) + 2;

    fmpq_init(b);
    arb_init(t);
    arb_init(f);

    bernoulli_fmpq_ui(b, n);
    arb_set_fmpq(x, b, wp);

    arb_const_pi(t, wp);
    arb_mul_2exp_si(t, t, 1);
    arb_pow_ui(t, t, n, wp);

    arb_fac_ui(f, n, wp);

    arb_div(t, t, f, wp);
    arb_mul(x, x, t, wp);
    arb_abs(x, x);
    arb_mul_2exp_si(x, x, -1);

    arb_clear(t);
    arb_clear(f);
    fmpq_clear(b);
}

