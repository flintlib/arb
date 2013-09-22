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

void fmprb_calc_newton_conv_factor(fmpr_t conv_factor,
    fmprb_calc_func_t func, void * param,
    const fmprb_t conv_region, long prec)
{
    fmprb_struct t[3];

    fmprb_init(t);
    fmprb_init(t + 1);
    fmprb_init(t + 2);

    func(t, conv_region, param, 3, prec);

    fmprb_div(t, t + 2, t + 1, prec);
    fmprb_mul_2exp_si(t, t, -1);

    fmprb_get_abs_ubound_fmpr(conv_factor, t, prec);

    fmprb_clear(t);
    fmprb_clear(t + 1);
    fmprb_clear(t + 2);
}

