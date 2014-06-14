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
#include "acb_poly.h"

void
arb_const_glaisher_eval(arb_t y, long prec)
{
    acb_struct z[2];
    acb_t s, a;
    long wp;

    acb_init(z + 0);
    acb_init(z + 1);
    acb_init(s);
    acb_init(a);

    wp = prec + 20;

    /* directly evaluating at s = -1 is slightly faster
       than evaluating at s = 2 */
    acb_set_si(s, -1);
    acb_one(a);
    _acb_poly_zeta_cpx_series(z, s, a, 0, 2, wp);

    arb_one(y);
    arb_div_ui(y, y, 12, wp);
    arb_sub(y, y, acb_realref(z + 1), wp);
    arb_exp(y, y, wp);

    acb_clear(z + 0);
    acb_clear(z + 1);
    acb_clear(s);
    acb_clear(a);
}

ARB_DEF_CACHED_CONSTANT(arb_const_glaisher, arb_const_glaisher_eval)

