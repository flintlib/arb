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

    Copyright (C) 2016 Pascal Molin

******************************************************************************/

#include "acb_dirichlet.h"
#include "acb_poly.h"

/* x = Pi / q * t^2 */
static void
acb_dirichlet_arb_theta_argt(arb_t x, ulong q, const arb_t t, slong prec)
{
    arb_const_pi(x, prec);
    arb_div_ui(x, x, q, prec);
    arb_mul(x, x, t, prec);
    arb_mul(x, x, t, prec);
    arb_neg(x, x);
    arb_exp(x, x, prec);
}

void
acb_dirichlet_chi_theta_arb(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, const arb_t t, slong prec)
{
    slong len;
    ulong * a;
    arb_t x;
    acb_dirichlet_powers_t z;

    len = acb_dirichlet_theta_length(G->q, t, prec);

    a = flint_malloc(len * sizeof(ulong));
    acb_dirichlet_ui_chi_vec(a, G, chi, len);
    acb_dirichlet_powers_init(z, chi->order.n, len, prec);

    arb_init(x);
    acb_dirichlet_arb_theta_argt(x, G->q, t, prec);
    acb_dirichlet_arb_theta_naive(res, x, chi->parity, a, z, len, prec);

    arb_clear(x);
    flint_free(a);
    acb_dirichlet_powers_clear(z);
}
