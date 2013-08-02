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

void
fmprb_const_glaisher_eval(fmprb_t y, long prec)
{
    fmprb_t t;
    fmprb_struct z[2];
    long wp;

    fmprb_init(t);
    fmprb_init(z + 0);
    fmprb_init(z + 1);

    wp = prec + 20;

    fmprb_set_ui(z + 0, 2);
    fmprb_one(z + 1);
    fmprb_one(t);

    _fmprb_poly_zeta_series(z, z, 2, t, 0, 2, wp);

    fmprb_const_pi(z, wp);
    fmprb_mul(z, z, z, wp);
    fmprb_mul_2exp_si(z, z, 1);
    fmprb_div(z, z + 1, z, wp);

    fmprb_const_euler(t, wp);
    fmprb_div_ui(t, t, 12, wp);
    fmprb_sub(t, t, z, wp);
    fmprb_exp(y, t, wp);

    fmprb_const_pi(t, wp);
    fmprb_mul_2exp_si(t, t, 1);
    fmprb_root(t, t, 12, wp);

    fmprb_mul(y, y, t, prec);

    fmprb_clear(t);
    fmprb_clear(z + 0);
    fmprb_clear(z + 1);
}

DEF_CACHED_CONSTANT(fmprb_const_glaisher, fmprb_const_glaisher_eval)

