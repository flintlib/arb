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

#include "fmprb_poly.h"

long _fmpr_mag(const fmpr_t c)
{
    long m = fmpz_bits(fmpr_manref(c)) + fmpz_get_si(fmpr_expref(c));

    return FLINT_MAX(m, 0);
}

void
_fmprb_poly_newton_refine_root(fmprb_t r, fmprb_srcptr poly, long len,
    const fmprb_t start,
    const fmprb_t convergence_interval,
    const fmpr_t convergence_factor,
    long eval_extra_prec,
    long prec)
{
    long precs[FLINT_BITS];
    long i, iters, wp, padding, start_prec;

    start_prec = fmprb_rel_accuracy_bits(start);

    padding = 5 + _fmpr_mag(convergence_factor);
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

    fmprb_set(r, start);

    for (i = iters - 1; i >= 0; i--)
    {
        wp = precs[i] + eval_extra_prec;

        if (!_fmprb_poly_newton_step(r, poly, len, r, convergence_interval,
            convergence_factor, wp))
        {
            printf("warning: newton_refine_root: improvement failed\n");
            break;
        }

    }
}

