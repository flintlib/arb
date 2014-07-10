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

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include "arb.h"
#include "acb.h"

int polylog_is_real(const acb_t s, const acb_t z);

void
arb_polylog(arb_t w, const arb_t s, const arb_t z, long prec)
{
    acb_t ss, zz;
    acb_init(ss);
    acb_init(zz);
    acb_set_arb(ss, s);
    acb_set_arb(zz, z);
    if (polylog_is_real(ss, zz))
    {
        acb_polylog(zz, ss, zz, prec);
        arb_set(w, acb_realref(zz));
    }
    else
    {
        arb_indeterminate(w);
    }
    acb_clear(ss);
    acb_clear(zz);
}

void
arb_polylog_si(arb_t w, long s, const arb_t z, long prec)
{
    arb_t t;
    arb_init(t);
    arb_set_si(t, s);
    arb_polylog(w, t, z, prec);
    arb_clear(t);
}
