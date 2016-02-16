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

#include "arb.h"
#include "acb_dirichlet.h"

void
arb_zeta_inv_ui_euler_product(arb_t z, ulong s, slong prec)
{
    const signed char chi[1] = {1};
    _acb_dirichlet_euler_product_real_ui(z, s, chi, 1, 1, prec);
}

void
arb_zeta_ui_euler_product(arb_t z, ulong s, slong prec)
{
    const signed char chi[1] = {1};
    _acb_dirichlet_euler_product_real_ui(z, s, chi, 1, 0, prec);
}

void
arb_zeta_ui_asymp(arb_t z, ulong s, slong prec)
{
    const signed char chi[1] = {1};
    _acb_dirichlet_euler_product_real_ui(z, s, chi, 1, 0, prec);
}

