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

#include "arb_poly.h"

static __inline__ long _arf_mag(const arf_t c)
{
    long m = arf_abs_bound_lt_2exp_si(c);
    return FLINT_MAX(m, 0);
}

void
_arb_poly_newton_refine_root(arb_t r, arb_srcptr poly, long len,
    const arb_t start,
    const arb_t convergence_interval,
    const arf_t convergence_factor,
    long eval_extra_prec,
    long prec)
{
    long precs[FLINT_BITS];
    long i, iters, wp, padding, start_prec;

    start_prec = arb_rel_accuracy_bits(start);

    padding = 5 + _arf_mag(convergence_factor);
    precs[0] = prec + padding;
    iters = 1;
    while ((iters < FLINT_BITS) && (precs[iters-1] + padding > 2*start_prec))
    {
        precs[iters] = (precs[iters-1] / 2) + padding;
        iters++;

        if (iters == FLINT_BITS)
        {
            printf("newton_refine_root: initial value too imprecise\n");
            abort();
        }
    }

    arb_set(r, start);

    for (i = iters - 1; i >= 0; i--)
    {
        wp = precs[i] + eval_extra_prec;

        if (!_arb_poly_newton_step(r, poly, len, r, convergence_interval,
            convergence_factor, wp))
        {
            printf("warning: newton_refine_root: improvement failed\n");
            break;
        }

    }
}

