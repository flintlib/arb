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

#include "fmprb_calc.h"

int fmprb_calc_refine_root_newton(fmprb_t r, fmprb_calc_func_t func,
    void * param, const fmprb_t start, const fmprb_t conv_region,
    const fmpr_t conv_factor, long eval_extra_prec, long prec)
{
    long precs[FLINT_BITS];
    long i, iters, wp, padding, start_prec;
    int result;

    start_prec = fmprb_rel_accuracy_bits(start);

    if (fmprb_calc_verbose)
        printf("newton initial accuracy: %ld\n", start_prec);

    padding = fmpr_abs_bound_lt_2exp_si(conv_factor);
    padding = FLINT_MIN(padding, prec) + 5;
    padding = FLINT_MAX(0, padding);

    precs[0] = prec + padding;
    iters = 1;
    while ((iters < FLINT_BITS) && (precs[iters-1] + padding > 2*start_prec))
    {
        precs[iters] = (precs[iters-1] / 2) + padding;
        iters++;

        if (iters == FLINT_BITS)
        {
            return FMPRB_CALC_IMPRECISE_INPUT;
        }
    }

    fmprb_set(r, start);

    for (i = iters - 1; i >= 0; i--)
    {
        wp = precs[i] + eval_extra_prec;

        if (fmprb_calc_verbose)
            printf("newton step: wp = %ld + %ld = %ld\n",
                precs[i], eval_extra_prec, wp);

        if ((result = fmprb_calc_newton_step(r, func, param,
            r, conv_region, conv_factor, wp)) != FMPRB_CALC_SUCCESS)
        {
            return result;
        }
    }

    return FMPRB_CALC_SUCCESS;
}

