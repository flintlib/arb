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

#include "arb.h"

void arb_zeta_inv_ui_euler_product(arb_t z, ulong s, long prec);

void
arb_bernoulli_ui_zeta(arb_t b, ulong n, long prec)
{
    long wp, piwp;

    arb_t t, u;

    if (n < 10 || n % 2 != 0)
        abort();

    wp = prec + 8;
    piwp = wp + 2*FLINT_BIT_COUNT(n);

    arb_init(t);
    arb_init(u);

    /* |B_n| = 2 * n! / (2*pi)^n * zeta(n) */
    arb_fac_ui(b, n, piwp);
    arb_const_pi(t, piwp);
    arb_mul_2exp_si(t, t, 1);
    arb_pow_ui(t, t, n, piwp);

    if (n > 0.7 * wp)
    {
        arb_zeta_ui_asymp(u, n, wp);
        arb_mul(b, b, u, wp);
    }
    else
    {
        arb_zeta_inv_ui_euler_product(u, n, wp);
        arb_mul(t, t, u, wp);
    }

    arb_div(b, b, t, prec);
    arb_mul_2exp_si(b, b, 1);

    if (n % 4 == 0)
        arb_neg(b, b);

    arb_clear(t);
    arb_clear(u);
}

