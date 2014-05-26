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

#include "arb.h"

void
arb_sinh(arb_t s, const arb_t x, long prec)
{
    if (arb_is_zero(x))
    {
        arb_zero(s);
    }
    else
    {
        arb_t t;
        long wp = prec + 4;

        arb_init(t);

        if (arf_cmpabs_2exp_si(arb_midref(x), -1) <= 0)
        {
            arb_mul_2exp_si(s, x, 1);
            arb_expm1(s, s, wp);
            arb_add_ui(t, s, 1, wp);
            arb_sqrt(t, t, wp);
            arb_div(s, s, t, prec);
        }
        else
        {
            arb_exp(s, x, wp);
            arb_inv(t, s, wp);
            arb_sub(s, s, t, prec);
        }

        arb_mul_2exp_si(s, s, -1);
        arb_clear(t);
    }
}

void
arb_cosh(arb_t c, const arb_t x, long prec)
{
    if (arb_is_zero(x))
    {
        arb_one(c);
    }
    else
    {
        arb_t t;
        long wp = prec + 4;

        arb_init(t);

        arb_exp(c, x, wp);
        arb_inv(t, c, wp);
        arb_add(c, c, t, prec);
        arb_mul_2exp_si(c, c, -1);

        arb_clear(t);
    }
}

void
arb_sinh_cosh(arb_t s, arb_t c, const arb_t x, long prec)
{
    if (arb_is_zero(x))
    {
        arb_zero(s);
        arb_one(c);
    }
    else
    {
        long wp = prec + 4;

        arb_t t, u;
        arb_init(t);
        arb_init(u);

        if (arf_cmpabs_2exp_si(arb_midref(x), -1) <= 0)
        {
            arb_mul_2exp_si(t, x, 1);
            arb_expm1(t, t, wp);
            arb_add_ui(u, t, 1, wp);
            arb_sqrt(u, u, wp);
            arb_inv(u, u, wp);
            arb_mul(s, t, u, prec);
            arb_add_ui(t, t, 2, wp);
            arb_mul(c, t, u, prec);
        }
        else
        {
            arb_exp(t, x, wp);
            arb_inv(u, t, wp);
            arb_sub(s, t, u, prec);
            arb_add(c, t, u, prec);
        }

        arb_mul_2exp_si(s, s, -1);
        arb_mul_2exp_si(c, c, -1);
        arb_clear(t);
        arb_clear(u);
    }
}

