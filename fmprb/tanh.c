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
fmprb_tanh(fmprb_t y, const fmprb_t x, long prec)
{
    fmprb_t t, u;

    fmprb_init(t);
    fmprb_init(u);

    fmprb_mul_2exp_si(t, x, 1);

    if (fmpr_sgn(fmprb_midref(x)) >= 0)
    {
        fmprb_neg(t, t);
        fmprb_expm1(t, t, prec + 4);
        fmprb_add_ui(y, t, 2, prec + 4);
        fmprb_div(y, t, y, prec);
        fmprb_neg(y, y);
    }
    else
    {
        fmprb_expm1(t, t, prec + 4);
        fmprb_add_ui(y, t, 2, prec + 4);
        fmprb_div(y, t, y, prec);
    }

    fmprb_clear(t);
    fmprb_clear(u);
}

