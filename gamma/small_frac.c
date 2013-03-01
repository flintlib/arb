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

#include "gamma.h"
#include "hypgeom.h"

void
gamma_const_1_3_eval(fmprb_t s, long prec)
{
    hypgeom_t series;
    fmprb_t t, u;

    fmprb_init(t);
    fmprb_init(u);

    hypgeom_init(series);

    fmpz_poly_set_str(series->A, "1  1");
    fmpz_poly_set_str(series->B, "1  1");
    fmpz_poly_set_str(series->P, "4  5 -46 108 -72");
    fmpz_poly_set_str(series->Q, "4  0 0 0 512000");

    prec += FLINT_CLOG2(prec);
    fmprb_hypgeom_infsum(s, t, series, prec, prec);

    fmprb_sqrt_ui(u, 10, prec);
    fmprb_mul(t, t, u, prec);

    fmprb_const_pi(u, prec);
    fmprb_pow_ui(u, u, 4, prec);
    fmprb_mul_ui(u, u, 12, prec);
    fmprb_mul(s, s, u, prec);

    fmprb_div(s, s, t, prec);
    fmprb_root(s, s, 2, prec);
    fmprb_root(s, s, 3, prec);

    hypgeom_clear(series);
    fmprb_clear(t);
    fmprb_clear(u);
}

DEF_CACHED_CONSTANT(gamma_const_1_3, gamma_const_1_3_eval)

void
gamma_const_1_4_eval(fmprb_t x, long prec)
{
    fmprb_t t, u;
    long wp = prec + 4 + 2 * FLINT_BIT_COUNT(prec);

    fmprb_init(t);
    fmprb_init(u);

    fmprb_one(t);
    fmprb_sqrt_ui(u, 2, wp);
    fmprb_agm(x, t, u, wp);

    fmprb_const_pi(t, wp);
    fmprb_mul_2exp_si(t, t, 1);
    fmprb_sqrt(u, t, wp);
    fmprb_mul(u, u, t, wp);

    fmprb_div(x, u, x, wp);
    fmprb_sqrt(x, x, wp);

    fmprb_clear(t);
    fmprb_clear(u);
}

DEF_CACHED_CONSTANT(gamma_const_1_4, gamma_const_1_4_eval)

void
gamma_small_frac(fmprb_t y, unsigned int p, unsigned int q, long prec)
{
    long wp = prec + 4;

    if (q == 1)
    {
        fmprb_one(y);
    }
    else if (q == 2)  /* p = 1 */
    {
        fmprb_const_sqrt_pi(y, prec);
    }
    else if (q == 3)
    {
        if (p == 1)
        {
            gamma_const_1_3(y, prec);
        }
        else  /* p = 2 */
        {
            fmprb_t t;
            fmprb_init(t);
            gamma_const_1_3(y, wp);
            fmprb_sqrt_ui(t, 3, wp);
            fmprb_mul(y, y, t, wp);
            fmprb_const_pi(t, wp);
            fmprb_div(y, t, y, prec);
            fmprb_mul_2exp_si(y, y, 1);
            fmprb_clear(t);
        }
    }
    else if (q == 4)
    {
        if (p == 1)
        {
            gamma_const_1_4(y, prec);
        }
        else  /* p = 3 */
        {
            fmprb_t t;
            fmprb_init(t);
            fmprb_sqrt_ui(y, 2, wp);
            fmprb_const_pi(t, wp);
            fmprb_mul(y, y, t, wp);
            gamma_const_1_4(t, wp);
            fmprb_div(y, y, t, prec);
            fmprb_clear(t);
        }
    }
    else if (q == 6)
    {
        fmprb_t t;
        fmprb_init(t);
        fmprb_const_pi(t, wp);
        fmprb_div_ui(t, t, 3, wp);
        fmprb_sqrt(t, t, wp);
        fmprb_set_ui(y, 2);
        fmprb_root(y, y, 3, wp);
        fmprb_mul(t, t, y, wp);
        gamma_const_1_3(y, wp);
        fmprb_mul(y, y, y, prec);

        if (p == 1)
        {
            fmprb_div(y, y, t, prec);
        }
        else  /* p = 5 */
        {
            fmprb_div(y, t, y, wp);
            fmprb_const_pi(t, wp);
            fmprb_mul(y, y, t, prec);
            fmprb_mul_2exp_si(y, y, 1);
        }

        fmprb_clear(t);
    }
    else
    {
        printf("small fraction not implemented!\n");
        abort();
    }
}

