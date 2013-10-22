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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "fmpcb.h"

void
fmpcb_root_exp(fmpcb_t r, const fmpcb_t a, long m, long index, long prec)
{
    if (m == 1)
        fmpcb_set_round(r, a, prec);
    else if (m == -1)
        fmpcb_inv(r, a, prec);
    else
    {
        fmpcb_log(r, a, prec);

        if (index != 0)
        {
            fmprb_t t;
            fmprb_init(t);
            fmprb_const_pi(t, prec);
            fmprb_mul_2exp_si(t, t, 1);
            fmprb_mul_si(t, t, index, prec);
            fmprb_add(fmpcb_imagref(r), fmpcb_imagref(r), t, prec);
            fmprb_clear(t);
        }

        fmpcb_div_si(r, r, m, prec);
        fmpcb_exp(r, r, prec);
    }
}

void
fmpcb_root_newton(fmpcb_t r, const fmpcb_t a, long m, long index, long prec)
{
    if (m == 1)
        fmpcb_set_round(r, a, prec);
    else if (m == -1)
        fmpcb_inv(r, a, prec);
    else
    {
        long startprec;
        fmpcb_t t;

        fmpcb_init(t);

        startprec = 100;
        fmpcb_root_exp(t, a, -FLINT_ABS(m), index, startprec);

        /* note: should check isolation first */
        fmpcb_invroot_newton(t, a, FLINT_ABS(m), t, startprec, prec);

        if (m < 0)
        {
            fmpcb_set(r, t);
        }
        else
        {
            if (fmpcb_is_one(a))
                fmpcb_conj(r, t);
            else if (m == 2)
                fmpcb_mul(r, t, a, prec);
            else
                fmpcb_inv(r, t, prec);
        }

        fmpcb_clear(t);
    }
}

void
fmpcb_root(fmpcb_t r, const fmpcb_t a, long m, long index, long prec)
{
    if (m == 1)
        fmpcb_set_round(r, a, prec);
    else if (m == -1)
        fmpcb_inv(r, a, prec);
    else if (m == 2)
        fmpcb_sqrt(r, a, prec);
    else if (m == -2)
        fmpcb_rsqrt(r, a, prec);
    else if (prec < 300 || !fmpcb_is_exact(a) || (m == LONG_MIN))
        fmpcb_root_exp(r, a, m, index, prec);
    else
        fmpcb_root_newton(r, a, m, index, prec);
}

