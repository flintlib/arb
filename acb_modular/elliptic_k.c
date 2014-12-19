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

#include "acb_modular.h"

void
acb_modular_elliptic_k(acb_t k, const acb_t m, long prec)
{
    acb_t t;
    acb_init(t);
    acb_sub_ui(t, m, 1, prec);
    acb_neg(t, t);
    acb_sqrt(t, t, prec);
    acb_agm1(k, t, prec);
    acb_const_pi(t, prec);
    acb_div(k, t, k, prec);
    acb_mul_2exp_si(k, k, -1);
    acb_clear(t);
}

