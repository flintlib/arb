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

#include "arith.h"
#include "zeta.h"

void
fmprb_zeta_ui_bernoulli(fmprb_t x, ulong n, long prec)
{
    fmpq_t b;
    fmprb_t t, f;
    long wp;

    if (n % 2)
        abort();

    wp = prec + FLINT_BIT_COUNT(n) + 2;

    fmpq_init(b);
    fmprb_init(t);
    fmprb_init(f);

    arith_bernoulli_number(b, n);
    fmprb_set_fmpq(x, b, wp);

    fmprb_const_pi(t, wp);
    fmprb_mul_2exp_si(t, t, 1);
    fmprb_pow_ui(t, t, n, wp);

    fmprb_fac_ui(f, n, wp);

    fmprb_div(t, t, f, wp);
    fmprb_mul(x, x, t, wp);
    fmprb_abs(x, x);
    fmprb_mul_2exp_si(x, x, -1);

    fmprb_clear(t);
    fmprb_clear(f);
    fmpq_clear(b);
}
