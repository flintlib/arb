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

#include "fmprb.h"
#include "fmprb_poly.h"
#include "arith.h"

void
_fmprb_cos_pi_fmpq_algebraic(fmprb_t c, ulong p, ulong q, long prec)
{
    /* handle simple angles using exact formulas */
    if (q <= 6)
    {
        if (p == 0)
        {
            fmprb_one(c);
        }
        else if (q == 2)  /* p/q must be 1/2 */
        {
            fmprb_zero(c);
        }
        else if (q == 3) /* p/q must be 1/3 */
        {
            fmprb_set_ui(c, 1);
            fmprb_mul_2exp_si(c, c, -1);
        }
        else if (q == 4)  /* p/q must be 1/4 */
        {
            fmprb_sqrt_ui(c, 2, prec);
            fmprb_mul_2exp_si(c, c, -1);
        }
        else if (q == 5) /* p/q must be 1/5 or 2/5 */
        {
            fmprb_sqrt_ui(c, 5, prec + 3);
            fmprb_add_si(c, c, (p == 1) ? 1 : -1, prec);
            fmprb_mul_2exp_si(c, c, -2);
        }
        else if (q == 6) /* p/q must be 1/6 */
        {
            fmprb_sqrt_ui(c, 3, prec);
            fmprb_mul_2exp_si(c, c, -1);
        }
    }
    /* reduce even denominator */
    else if (q % 2 == 0)
    {
        long extra = 2 * FLINT_BIT_COUNT(q) + 2;

        if (4 * p <= q)
        {
            _fmprb_cos_pi_fmpq_algebraic(c, p, q / 2, prec + extra);
            fmprb_add_ui(c, c, 1, prec + extra);
        }
        else
        {
            _fmprb_cos_pi_fmpq_algebraic(c, q / 2 - p, q / 2, prec + extra);
            fmprb_sub_ui(c, c, 1, prec + extra);
            fmprb_neg(c, c);
        }

        fmprb_mul_2exp_si(c, c, -1);
        fmprb_sqrt(c, c, prec);
    }
    else
    {
        /* compute root of the minimal polynomial */
        long start_prec, eval_extra_prec;
        fmpz_poly_t poly;
        fmprb_poly_t fpoly;
        fmpr_t interval_bound;
        fmprb_t interval;

        fmpr_init(interval_bound);
        fmprb_init(interval);
        fmpz_poly_init(poly);
        fmprb_poly_init(fpoly);

        if (p % 2 == 0)
            arith_cos_minpoly(poly, q);
        else
            arith_cos_minpoly(poly, 2 * q);

        eval_extra_prec = fmpz_poly_max_bits(poly);
        eval_extra_prec = FLINT_ABS(eval_extra_prec);
        fmprb_poly_set_fmpz_poly(fpoly, poly, FMPR_PREC_EXACT);

        /* todo: smallify for accuracy */
        start_prec = 100 + eval_extra_prec;
        fmprb_const_pi(c, start_prec);
        fmprb_mul_ui(c, c, p, start_prec);
        fmprb_div_ui(c, c, q, start_prec);
        fmprb_cos(c, c, start_prec);

        if (100 + eval_extra_prec - 10 < prec)
        {
            fmprb_set(interval, c);
            fmpr_mul_2exp_si(fmprb_radref(interval), fmprb_radref(interval), 1);
            _fmprb_poly_newton_convergence_factor(interval_bound,
                fpoly->coeffs, fpoly->length, interval, start_prec);
            _fmprb_poly_newton_refine_root(c, fpoly->coeffs, fpoly->length,
                c, interval, interval_bound, eval_extra_prec, prec);
        }

        fmpz_poly_clear(poly);
        fmprb_poly_clear(fpoly);
        fmpr_clear(interval_bound);
        fmprb_clear(interval);
    }
}

void
_fmprb_sin_pi_fmpq_algebraic(fmprb_t s, ulong p, ulong q, long prec)
{
    if (q % 2 == 0)
    {
        p = q / 2 - p;

        while ((p % 2 == 0) && (q % 2 == 0))
        {
            p /= 2;
            q /= 2;
        }

        _fmprb_cos_pi_fmpq_algebraic(s, p, q, prec);
    }
    else
    {
        _fmprb_cos_pi_fmpq_algebraic(s, q - 2 * p, 2 * q, prec);
    }
}

void
_fmprb_sin_cos_pi_fmpq_algebraic(fmprb_t s, fmprb_t c, ulong p, ulong q, long prec)
{
    long wp;

    if (q <= 6)
    {
        if (p == 0)
        {
            fmprb_one(c);
            fmprb_zero(s);
            return;
        }
        else if (q == 2)  /* p/q must be 1/2 */
        {
            fmprb_zero(c);
            fmprb_one(s);
            return;
        }
        else if (q == 4)  /* p/q must be 1/4 */
        {
            fmprb_sqrt_ui(c, 2, prec);
            fmprb_mul_2exp_si(c, c, -1);
            fmprb_set(s, c);
            return;
        }
    }

    wp = prec + 3;

    /* prefer the formula with less cancellation */
    if (p <= q / 4)
    {
        _fmprb_sin_pi_fmpq_algebraic(s, p, q, wp);
        fmprb_mul(c, s, s, wp);
        fmprb_sub_ui(c, c, 1, wp);
        fmprb_neg(c, c);
        fmprb_sqrt(c, c, prec);
    }
    else
    {
        _fmprb_cos_pi_fmpq_algebraic(c, p, q, wp);
        fmprb_mul(s, c, c, wp);
        fmprb_sub_ui(s, s, 1, wp);
        fmprb_neg(s, s);
        fmprb_sqrt(s, s, prec);
    }
}

