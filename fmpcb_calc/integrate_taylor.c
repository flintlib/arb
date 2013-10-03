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

#include "fmpcb_calc.h"
#include "math.h"

int
fmpcb_calc_integrate_taylor(fmpcb_t res,
    fmpcb_calc_func_t func, void * param,
    const fmpcb_t a, const fmpcb_t b,
    const fmpr_t inner_radius,
    const fmpr_t outer_radius,
    long accuracy_goal, long prec)
{
    long num_steps, step, N, bp;
    int result;

    fmpcb_t delta, m, x, y1, y2, sum;
    fmpcb_ptr taylor_poly;
    fmpr_t err;

    fmpcb_init(delta);
    fmpcb_init(m);
    fmpcb_init(x);
    fmpcb_init(y1);
    fmpcb_init(y2);
    fmpcb_init(sum);
    fmpr_init(err);

    fmpcb_sub(delta, b, a, prec);

    /* precision used for bounds calculations */
    bp = FMPRB_RAD_PREC;

    /* compute the number of steps */
    {
        fmpr_t t;
        fmpr_init(t);
        fmpcb_get_abs_ubound_fmpr(t, delta, bp);
        fmpr_div(t, t, inner_radius, bp, FMPR_RND_UP);
        fmpr_mul_2exp_si(t, t, -1);
        num_steps = (long) (fmpr_get_d(t, FMPR_RND_UP) + 1.0);
        /* make sure it's not something absurd */
        num_steps = FLINT_MIN(num_steps, 10 * prec);
        num_steps = FLINT_MAX(num_steps, 1);
        fmpr_clear(t);
    }

    result = FMPRB_CALC_SUCCESS;

    fmpcb_zero(sum);

    for (step = 0; step < num_steps; step++)
    {
        /* midpoint of subinterval */
        fmpcb_mul_ui(m, delta, 2 * step + 1, prec);
        fmpcb_div_ui(m, m, 2 * num_steps, prec);
        fmpcb_add(m, m, a, prec);

        if (fmprb_calc_verbose)
        {
            printf("integration point %ld/%ld: ", 2 * step + 1, 2 * num_steps);
            fmpcb_printd(m, 15); printf("\n");
        }

        /* evaluate at +/- x */
        /* TODO: exactify m, and include error in x? */
        fmpcb_div_ui(x, delta, 2 * num_steps, prec);

        /* compute bounds and number of terms to use */
        {
            fmprb_t cbound, xbound, rbound;
            fmpr_t C, D, R, X, T;
            double DD, TT, NN;

            fmprb_init(cbound);
            fmprb_init(xbound);
            fmprb_init(rbound);
            fmpr_init(C);
            fmpr_init(D);
            fmpr_init(R);
            fmpr_init(X);
            fmpr_init(T);

            /* R is the outer radius */
            fmpr_set(R, outer_radius);

            /* X = upper bound for |x| */
            fmpcb_get_abs_ubound_fmpr(X, x, bp);
            fmprb_set_fmpr(xbound, X);

            /* Compute C(m,R). Important subtlety: due to rounding when
               computing m, we will in general be farther than R away from
               the integration path. But since fmpcb_calc_cauchy_bound
               actually integrates over the area traced by a complex
               interval, it will catch any extra singularities (giving
               an infinite bound). */
            fmprb_set_fmpr(rbound, outer_radius);
            fmpcb_calc_cauchy_bound(cbound, func, param, m, rbound, 8, bp);
            fmpr_add(C, fmprb_midref(cbound), fmprb_radref(cbound), bp, FMPR_RND_UP);

            /* Sanity check: we need C < inf and R > X */
            if (fmpr_is_finite(C) && fmpr_cmp(R, X) > 0)
            {
                /* Compute upper bound for D = C * R * X / (R - X) */
                fmpr_mul(D, C, R, bp, FMPR_RND_UP);
                fmpr_mul(D, D, X, bp, FMPR_RND_UP);
                fmpr_sub(T, R, X, bp, FMPR_RND_DOWN);
                fmpr_div(D, D, T, bp, FMPR_RND_UP);

                /* Compute upper bound for T = (X / R) */
                fmpr_div(T, X, R, bp, FMPR_RND_UP);

                /* Choose N */
                /* TODO: use fmpr arithmetic to avoid overflow */
                /* TODO: use relative accuracy (look at |f(m)|?) */
                DD = fmpr_get_d(D, FMPR_RND_UP);
                TT = fmpr_get_d(T, FMPR_RND_UP);
                NN = -(accuracy_goal * 0.69314718055994530942 + log(DD)) / log(TT);
                N = NN + 0.5;
                N = FLINT_MIN(N, 100 * prec);
                N = FLINT_MAX(N, 1);

                /* Tail bound: D / (N + 1) * T^N */
                fmpr_pow_sloppy_ui(T, T, N, bp, FMPR_RND_UP);
                fmpr_mul(D, D, T, bp, FMPR_RND_UP);
                fmpr_div_ui(err, D, N + 1, bp, FMPR_RND_UP);
            }
            else
            {
                N = 1;
                fmpr_pos_inf(err);
                result = FMPRB_CALC_NO_CONVERGENCE;
            }

            if (fmprb_calc_verbose)
            {
                printf("N = %ld; bound: ", N); fmpr_printd(err, 15); printf("\n");
                printf("R: "); fmpr_printd(R, 15); printf("\n");
                printf("C: "); fmpr_printd(C, 15); printf("\n");
                printf("X: "); fmpr_printd(X, 15); printf("\n");
            }

            fmprb_clear(cbound);
            fmprb_clear(xbound);
            fmprb_clear(rbound);
            fmpr_clear(C);
            fmpr_clear(D);
            fmpr_clear(R);
            fmpr_clear(X);
            fmpr_clear(T);
        }

        /* evaluate Taylor polynomial */
        taylor_poly = _fmpcb_vec_init(N + 1);
        func(taylor_poly, m, param, N, prec);
        _fmpcb_poly_integral(taylor_poly, taylor_poly, N + 1, prec);
        _fmpcb_poly_evaluate(y2, taylor_poly, N + 1, x, prec);
        fmpcb_neg(x, x);
        _fmpcb_poly_evaluate(y1, taylor_poly, N + 1, x, prec);
        fmpcb_neg(x, x);

        /* add truncation error */
        fmprb_add_error_fmpr(fmpcb_realref(y1), err);
        fmprb_add_error_fmpr(fmpcb_imagref(y1), err);
        fmprb_add_error_fmpr(fmpcb_realref(y2), err);
        fmprb_add_error_fmpr(fmpcb_imagref(y2), err);

        fmpcb_add(sum, sum, y2, prec);
        fmpcb_sub(sum, sum, y1, prec);

        if (fmprb_calc_verbose)
        {
            printf("values:  ");
            fmpcb_printd(y1, 15); printf("  ");
            fmpcb_printd(y2, 15); printf("\n");
        }

        _fmpcb_vec_clear(taylor_poly, N + 1);

        if (result == FMPRB_CALC_NO_CONVERGENCE)
            break;
    }

    fmpcb_set(res, sum);

    fmpcb_clear(delta);
    fmpcb_clear(m);
    fmpcb_clear(x);
    fmpcb_clear(y1);
    fmpcb_clear(y2);
    fmpcb_clear(sum);
    fmpr_clear(err);

    return result;
}

