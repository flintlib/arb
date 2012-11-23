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

/*
Sets r = exp((-1)^inverse * (1/m) * (log(a) + 2 pi i k)). As k goes from 0 to m-1, this
expression gives all the mth (inverse) roots of the complex number a, starting with
the principal mth root.

Algorithm: if the precision is high enough and |m| is small enough, we
use Newton iteration to compute an inverse |m|-th root by solving
f(z) = (1/z)^m - a = 0, by calling fmpcb_invroot_newton.
The initial value is obtained by evaluating the exponential function.
TODO: we should check that the initial value isolates the right root
(if the error of a is large).
*/

void
fmpcb_root_via_exp(fmpcb_t r, const fmpcb_t a, ulong k, ulong m, int inverse, long prec)
{
    fmpcb_log(r, a, prec);

    if (k != 0)
    {
        fmprb_t t;
        fmprb_init(t);
        fmprb_const_pi(t, prec);
        fmprb_mul_2exp_si(t, t, 1);
        fmprb_mul_ui(t, t, k, prec);
        fmprb_add(fmpcb_imagref(r), fmpcb_imagref(r), t, prec);
        fmprb_clear(t);
    }

    fmpcb_div_ui(r, r, m, prec);
    if (inverse)
        fmpcb_neg(r, r);

    fmpcb_exp(r, r, prec);
}

/* fixme: m+2 overflows in Newton */



void
fmpcb_root(fmpcb_t r, const fmpcb_t a, ulong k, ulong m, int inverse, long prec)
{
    if (m == 1)
    {
        if (inverse)
            fmpcb_inv(r, a, prec);
        else
            fmpcb_set_round(r, a, prec);
    }
    else if (prec < 300)
    {
        fmpcb_root_via_exp(r, a, k, m, inverse, prec);
    }
    else
    {
        long startprec;

        fmpcb_t t;
        fmpcb_init(t);

        startprec = 100;
        fmpcb_root_via_exp(t, a, k, m, 1, startprec);

        /* note: should check isolation first */
        if (1)
        {
            fmpcb_invroot_newton(t, a, m, t, startprec, prec);
        }

        if (inverse)
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

