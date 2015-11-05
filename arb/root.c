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

    Copyright (C) 2015 Fredrik Johansson

******************************************************************************/

#include "arb.h"

static void
arb_root_arf(arb_t z, const arf_t x, ulong k, slong prec)
{
    int inexact = arf_root(arb_midref(z), x, k, prec, ARB_RND);

    if (inexact)
        arf_mag_set_ulp(arb_radref(z), arb_midref(z), prec);
    else
        mag_zero(arb_radref(z));
}

void
arb_root_algebraic(arb_t res, const arb_t x, ulong k, slong prec)
{
    mag_t r, msubr, m1k, t;

    if (arb_is_exact(x))
    {
        arb_root_arf(res, arb_midref(x), k, prec);
        return;
    }

    if (!arb_is_nonnegative(x))
    {
        arb_indeterminate(res);
        return;
    }

    mag_init(r);
    mag_init(msubr);
    mag_init(m1k);
    mag_init(t);

    /* x = [m-r, m+r] */
    mag_set(r, arb_radref(x));
    /* m - r */
    arb_get_mag_lower(msubr, x);

    /* m^(1/k) */
    arb_root_arf(res, arb_midref(x), k, prec);

    /* bound for m^(1/k) */
    arb_get_mag(m1k, res);

    /* C = min(1, log(1+r/(m-r))/k) */
    mag_div(t, r, msubr);
    mag_log1p(t, t);
    mag_div_ui(t, t, k);
    if (mag_cmp_2exp_si(t, 0) > 0)
        mag_one(t);

    /* C m^(1/k) */
    mag_mul(t, m1k, t);
    mag_add(arb_radref(res), arb_radref(res), t);

    mag_clear(r);
    mag_clear(msubr);
    mag_clear(m1k);
    mag_clear(t);
}

void
arb_root_exp(arb_t res, const arb_t x, ulong k, slong prec)
{
    arb_log(res, x, prec + 4);
    arb_div_ui(res, res, k, prec + 4);
    arb_exp(res, res, prec);
}

void
arb_root(arb_t res, const arb_t x, ulong k, slong prec)
{
    if (k == 0)
    {
        arb_indeterminate(res);
    }
    else if (k == 1)
    {
        arb_set_round(res, x, prec);
    }
    else if (k == 2)
    {
        arb_sqrt(res, x, prec);
    }
    else if (k == 4)
    {
        arb_sqrt(res, x, prec + 2);
        arb_sqrt(res, res, prec);
    }
    else
    {
        if (k > 50 || prec < (1L << ((k / 8) + 8)))
            arb_root_exp(res, x, k, prec);
        else
            arb_root_algebraic(res, x, k, prec);
    }
}

