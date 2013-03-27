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
#include "fmpcb.h"
#include "zeta.h"

void
fmprb_zeta(fmprb_t y, const fmprb_t s, long prec)
{
    fmpcb_t t;
    fmpcb_init(t);
    fmpcb_set_fmprb(t, s),
    fmpcb_zeta(t, t, prec);
    fmprb_set(y, fmpcb_realref(t));
    fmpcb_clear(t);
}

void
fmprb_zeta_ui(fmprb_t b, ulong n, long prec)
{
    zeta_ui(b, n, prec);
}

