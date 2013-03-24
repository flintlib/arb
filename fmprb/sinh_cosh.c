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

void
fmprb_sinh(fmprb_t s, const fmprb_t x, long prec)
{
    if (fmprb_is_zero(x))
    {
        fmprb_zero(s);
    }
    else
    {
        fmprb_t t;
        long wp = prec + 4;

        fmprb_init(t);

        if (fmpr_cmpabs_2exp_si(fmprb_midref(x), -1) <= 0)
        {
            fmprb_mul_2exp_si(s, x, 1);
            fmprb_expm1(s, s, wp);
            fmprb_add_ui(t, s, 1, wp);
            fmprb_sqrt(t, t, wp);
            fmprb_div(s, s, t, prec);
        }
        else
        {
            fmprb_exp(s, x, wp);
            fmprb_ui_div(t, 1, s, wp);
            fmprb_sub(s, s, t, prec);
        }

        fmprb_mul_2exp_si(s, s, -1);
        fmprb_clear(t);
    }
}

void
fmprb_cosh(fmprb_t c, const fmprb_t x, long prec)
{
    if (fmprb_is_zero(x))
    {
        fmprb_one(c);
    }
    else
    {
        fmprb_t t;
        long wp = prec + 4;

        fmprb_init(t);

        fmprb_exp(c, x, wp);
        fmprb_ui_div(t, 1, c, wp);
        fmprb_add(c, c, t, prec);
        fmprb_mul_2exp_si(c, c, -1);

        fmprb_clear(t);
    }
}

void
fmprb_sinh_cosh(fmprb_t s, fmprb_t c, const fmprb_t x, long prec)
{
    if (fmprb_is_zero(x))
    {
        fmprb_zero(s);
        fmprb_one(c);
    }
    else
    {
        long wp = prec + 4;

        fmprb_t t, u;
        fmprb_init(t);
        fmprb_init(u);

        if (fmpr_cmpabs_2exp_si(fmprb_midref(x), -1) <= 0)
        {
            fmprb_mul_2exp_si(t, x, 1);
            fmprb_expm1(t, t, wp);
            fmprb_add_ui(u, t, 1, wp);
            fmprb_sqrt(u, u, wp);
            fmprb_ui_div(u, 1, u, wp);
            fmprb_mul(s, t, u, prec);
            fmprb_add_ui(t, t, 2, wp);
            fmprb_mul(c, t, u, prec);
        }
        else
        {
            fmprb_exp(t, x, wp);
            fmprb_ui_div(u, 1, t, wp);
            fmprb_sub(s, t, u, prec);
            fmprb_add(c, t, u, prec);
        }

        fmprb_mul_2exp_si(s, s, -1);
        fmprb_mul_2exp_si(c, c, -1);
        fmprb_clear(t);
        fmprb_clear(u);
    }
}

