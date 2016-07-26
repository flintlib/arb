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

/* x(t) = exp(-Pi / q * t^2) */
static void
acb_dirichlet_arb_theta_xt(arb_t xt, ulong q, const arb_t t, slong prec)
{
    arb_const_pi(xt, prec);
    arb_div_ui(xt, xt, q, prec);
    arb_mul(xt, xt, t, prec);
    arb_mul(xt, xt, t, prec);
    arb_neg(xt, xt);
    arb_exp(xt, xt, prec);
}

void
acb_dirichlet_chi_theta_arb_smallorder(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, const arb_t xt, slong len, slong prec)
{
    ulong * a;
    acb_dirichlet_powers_t z;
    
    a = flint_malloc(len * sizeof(ulong));
    acb_dirichlet_ui_chi_vec(a, G, chi, len);

    _acb_dirichlet_powers_init(z, chi->order.n, 0, 0, prec);
    acb_dirichlet_arb_theta_smallorder(res, xt, chi->parity, a, z, len, prec);
    acb_dirichlet_powers_clear(z);

    flint_free(a);
}

void
acb_dirichlet_chi_theta_arb_series(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, const arb_t xt, slong len, slong prec)
{
    acb_ptr a;
    a = _acb_vec_init(len);
    acb_dirichlet_chi_vec(a, G, chi, len, prec);
    if (chi->parity)
    {
        slong k;
        for (k = 2; k < len; k++)
            acb_mul_si(a + k, a + k, k, prec);
    }
    acb_dirichlet_qseries_eval_arb(res, a, xt, len, prec);
    _acb_vec_clear(a, len);
}

void
acb_dirichlet_vec_mellin_series(acb_ptr res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong len, const arb_t t, slong n, slong prec)
{
    slong k;
    arb_t tk, xt, stk, st;
    acb_ptr a;
    a = _acb_vec_init(len);
    acb_dirichlet_chi_vec(a, G, chi, len, prec);
    if (chi->parity)
    {
        for (k = 2; k < len; k++)
            acb_mul_si(a + k, a + k, k, prec);
    }
    arb_init(tk);
    arb_init(xt);
    arb_init(st);
    arb_init(stk);

    arb_sqrt(st, t, prec);
    arb_one(tk);
    arb_one(stk);
    for (k = 0; k < n; k++)
    {
        acb_dirichlet_arb_theta_xt(xt, G->q, tk, prec);
        /* TODO: reduce len */
        acb_dirichlet_qseries_eval_arb(res + k, a, xt, len, prec);
        acb_mul_arb(res + k, res + k, stk, prec);
        arb_mul(tk, tk, t, prec);
        arb_mul(stk, stk, st, prec);
    }
    arb_clear(xt);
    arb_clear(tk);
    arb_clear(stk);
    arb_clear(st);
    _acb_vec_clear(a, len);
}

void
acb_dirichlet_chi_theta_arb_naive(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, const arb_t xt, slong len, slong prec)
{
    ulong * a;
    acb_dirichlet_powers_t z;

    a = flint_malloc(len * sizeof(ulong));
    acb_dirichlet_ui_chi_vec(a, G, chi, len);

    acb_dirichlet_powers_init(z, chi->order.n, len, prec);

    acb_dirichlet_arb_theta_naive(res, xt, chi->parity, a, z, len, prec);

    acb_dirichlet_powers_clear(z);
    flint_free(a);
}


void
acb_dirichlet_chi_theta_arb(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, const arb_t t, slong prec)
{
    slong len;
    arb_t xt;

    len = acb_dirichlet_theta_length(G->q, t, prec);

    arb_init(xt);
    acb_dirichlet_arb_theta_xt(xt, G->q, t, prec);

    /* TODO: tune this limit */
    if (chi->order.n < 30)
        acb_dirichlet_chi_theta_arb_smallorder(res, G, chi, xt, len, prec);
    else
        acb_dirichlet_chi_theta_arb_naive(res, G, chi, xt, len, prec);

    arb_clear(xt);
}
