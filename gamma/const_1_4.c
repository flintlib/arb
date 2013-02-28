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

