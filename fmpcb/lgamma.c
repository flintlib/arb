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

#include "fmpcb.h"
#include "gamma.h"

/* corrects branch cut of sum_{k=0}^{r-1} log(z+k), given the
   logarithm of the product */
void
_fmpcb_log_rising_correct_branch(fmpcb_t t,
        const fmpcb_t t_wrong, const fmpcb_t z, ulong r, long prec)
{
    fmpcb_t f;
    fmprb_t pi, u, v;
    fmpz_t pi_mult;
    long i, argprec;

    fmpcb_init(f);

    fmprb_init(u);
    fmprb_init(pi);
    fmprb_init(v);

    fmpz_init(pi_mult);

    argprec = FLINT_MIN(prec, 40);

    fmprb_zero(u);
    for (i = 0; i < r; i++)
    {
        fmpcb_add_ui(f, z, i, argprec);
        fmpcb_arg(v, f, argprec);
        fmprb_add(u, u, v, argprec);
    }

    if (argprec == prec)
    {
        fmprb_set(fmpcb_imagref(t), u);
    }
    else
    {
        fmprb_sub(v, u, fmpcb_imagref(t), argprec);
        fmprb_const_pi(pi, argprec);
        fmprb_div(v, v, pi, argprec);

        if (fmprb_get_unique_fmpz(pi_mult, v))
        {
            fmprb_const_pi(v, prec);
            fmprb_mul_fmpz(v, v, pi_mult, prec);
            fmprb_add(fmpcb_imagref(t), fmpcb_imagref(t), v, prec);
        }
        else
        {
            fmprb_zero(u);
            for (i = 0; i < r; i++)
            {
                fmpcb_add_ui(f, z, i, prec);
                fmpcb_arg(v, f, prec);
                fmprb_add(u, u, v, prec);
            }
            fmprb_set(fmpcb_imagref(t), u);
        }
    }

    fmpcb_clear(f);

    fmprb_clear(u);
    fmprb_clear(v);
    fmprb_clear(pi);

    fmpz_clear(pi_mult);
}

void
fmpcb_lgamma(fmpcb_t y, const fmpcb_t x, long prec)
{
    int reflect;
    long r, n, wp;
    fmpcb_t t, u;

    wp = prec + FLINT_BIT_COUNT(prec);

    gamma_stirling_choose_param_fmpcb(&reflect, &r, &n, x, 0, 0, wp);

    /* log(gamma(x)) = log(gamma(x+r)) - log(rf(x,r)) */
    fmpcb_init(t);
    fmpcb_init(u);

    fmpcb_add_ui(t, x, r, wp);
    gamma_stirling_eval_fmpcb(u, t, n, 0, wp);

    gamma_rising_fmpcb_ui_bsplit(t, x, r, prec);
    fmpcb_log(t, t, prec);

    _fmpcb_log_rising_correct_branch(t, t, x, r, wp);

    fmpcb_sub(y, u, t, prec);

    fmpcb_clear(t);
    fmpcb_clear(u);
}

