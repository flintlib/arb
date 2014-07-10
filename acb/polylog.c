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

#include "acb.h"
#include "acb_poly.h"

void
acb_polylog(acb_t w, const acb_t s, const acb_t z, long prec)
{
    acb_t t;
    acb_init(t);
    _acb_poly_polylog_cpx(t, s, z, 1, prec);
    acb_swap(w, t);
    acb_clear(t);
}

void
acb_polylog_si(acb_t w, long s, const acb_t z, long prec)
{
    acb_t t;
    acb_init(t);
    acb_set_si(t, s);
    acb_polylog(w, t, z, prec);
    acb_clear(t);
}
