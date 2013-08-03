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
#include "zeta.h"

void
fmprb_const_glaisher_eval(fmprb_t y, long prec)
{
    fmpcb_struct z[2];
    fmpcb_t s, a;
    long wp;

    fmpcb_init(z + 0);
    fmpcb_init(z + 1);
    fmpcb_init(s);
    fmpcb_init(a);

    wp = prec + 20;

    /* directly evaluating at s = -1 is slightly faster
       than evaluating at s = 2 */
    fmpcb_set_si(s, -1);
    fmpcb_one(a);
    zeta_series(z, s, a, 0, 2, wp);

    fmprb_one(y);
    fmprb_div_ui(y, y, 12, wp);
    fmprb_sub(y, y, fmpcb_realref(z + 1), wp);
    fmprb_exp(y, y, wp);

    fmpcb_clear(z + 0);
    fmpcb_clear(z + 1);
    fmpcb_clear(s);
    fmpcb_clear(a);
}

DEF_CACHED_CONSTANT(fmprb_const_glaisher, fmprb_const_glaisher_eval)

