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

/* TODO: refactor/combine some of this code with isolate_roots.c */

/* 0 means that it *could* be zero; otherwise +/- 1 */
static __inline__ int
_fmprb_sign(const fmprb_t t)
{
    if (fmprb_is_positive(t))
        return 1;
    else if (fmprb_is_negative(t))
        return -1;
    else
        return 0;
}

static int
partition(fmprb_t L, fmprb_t R,
    fmprb_calc_func_t func, void * param, const fmprb_t block, long prec)
{
    fmprb_t t, m;
    int msign;

    fmprb_init(t);
    fmprb_init(m);

    /* TODO: try other points */
    fmprb_set_fmpr(m, fmprb_midref(block));
    func(t, m, param, 1, prec);
    msign = _fmprb_sign(t);

    fmpr_mul_2exp_si(fmprb_radref(L), fmprb_radref(block), -1);
    fmpr_set(fmprb_radref(R), fmprb_radref(L));

    /* XXX: deal with huge shifts */
    fmpr_sub(fmprb_midref(L), fmprb_midref(block), fmprb_radref(L), FMPR_PREC_EXACT, FMPR_RND_DOWN);
    fmpr_add(fmprb_midref(R), fmprb_midref(block), fmprb_radref(R), FMPR_PREC_EXACT, FMPR_RND_DOWN);

    fmprb_clear(t);
    fmprb_clear(m);

    return msign;
}

int fmprb_calc_refine_root_bisect(fmprb_t r, fmprb_calc_func_t func,
    void * param, const fmprb_t start, long iter, long prec)
{
    int asign, bsign, msign, result;
    long i;
    fmprb_t t, u;

    fmprb_init(t);
    fmprb_init(u);

    /* XXX: deal with huge shifts */
    fmpr_sub(fmprb_midref(t), fmprb_midref(start), fmprb_radref(start),
        FMPR_PREC_EXACT, FMPR_RND_DOWN);
    func(u, t, param, 1, prec);
    asign = _fmprb_sign(u);

    fmpr_add(fmprb_midref(t), fmprb_midref(start), fmprb_radref(start),
        FMPR_PREC_EXACT, FMPR_RND_DOWN);
    func(u, t, param, 1, prec);
    bsign = _fmprb_sign(u);

    /* must have proper sign changes */
    if (asign == 0 || bsign == 0 || asign == bsign)
    {
        result = FMPRB_CALC_IMPRECISE_INPUT;
    }
    else
    {
        fmprb_set(r, start);

        result = FMPRB_CALC_SUCCESS;

        for (i = 0; i < iter; i++)
        {
            msign = partition(t, u, func, param, r, prec);

            /* the algorithm fails if the value at the midpoint cannot
               be distinguished from zero */
            if (msign == 0)
            {
                result = FMPRB_CALC_NO_CONVERGENCE;
                break;
            }

            if (msign == asign)
                fmprb_swap(r, u);
            else
                fmprb_swap(r, t);
        }
    }

    fmprb_clear(t);
    fmprb_clear(u);

    return result;
}

