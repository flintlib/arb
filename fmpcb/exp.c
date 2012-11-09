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
fmpcb_exp(fmpcb_t r, const fmpcb_t z, long prec)
{
#define a fmpcb_realref(z)
#define b fmpcb_imagref(z)

    fmprb_t t, u, v;

    fmprb_init(t);
    fmprb_init(u);
    fmprb_init(v);

    fmprb_exp(t, a, prec);
    fmprb_sin_cos(u, v, b, prec);

    fmprb_mul(fmpcb_realref(r), t, v, prec);
    fmprb_mul(fmpcb_imagref(r), t, u, prec);

    fmprb_clear(t);
    fmprb_clear(u);
    fmprb_clear(v);
}

