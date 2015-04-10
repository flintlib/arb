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

    Copyright (C) 2015 Fredrik Johansson

******************************************************************************/

#include "acb_hypgeom.h"

static void
_acb_hypgeom_li(acb_t res, const acb_t z, long prec)
{
    if (acb_is_zero(z))
    {
        acb_zero(res);
    }
    else
    {
        acb_log(res, z, prec);
        acb_hypgeom_ei(res, res, prec);
    }
}

void
_acb_hypgeom_const_li2_eval(arb_t s, long prec)
{
    acb_t t;
    acb_init(t);
    acb_set_ui(t, 2);
    _acb_hypgeom_li(t, t, prec);
    arb_set(s, acb_realref(t));
    acb_clear(t);
}

ARB_DEF_CACHED_CONSTANT(_acb_hypgeom_const_li2, _acb_hypgeom_const_li2_eval)

static void
_acb_hypgeom_li_offset(acb_t res, const acb_t z, long prec)
{
    if (acb_is_int(z) && arf_cmp_2exp_si(arb_midref(acb_realref(z)), 1) == 0)
    {
        acb_zero(res);
    }
    else
    {
        arb_t t;
        arb_init(t);
        _acb_hypgeom_const_li2(t, prec);
        _acb_hypgeom_li(res, z, prec);
        arb_sub(acb_realref(res), acb_realref(res), t, prec);
        arb_clear(t);
    }
}

void
acb_hypgeom_li(acb_t res, const acb_t z, int offset, long prec)
{
    if (offset)
        _acb_hypgeom_li_offset(res, z, prec);
    else
        _acb_hypgeom_li(res, z, prec);
}

