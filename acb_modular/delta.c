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
acb_modular_delta(acb_t z, const acb_t tau, long prec)
{
    psl2z_t g;
    arf_t one_minus_eps;
    acb_t tau_prime, t1, t2, t3, t4, w, q;

    psl2z_init(g);
    arf_init(one_minus_eps);
    acb_init(tau_prime);
    acb_init(t1);
    acb_init(t2);
    acb_init(t3);
    acb_init(t4);
    acb_init(w);
    acb_init(q);

    arf_set_ui_2exp_si(one_minus_eps, 63, -6);
    acb_modular_fundamental_domain_approx(tau_prime, g, tau,
        one_minus_eps, prec);

    acb_one(w);
    acb_exp_pi_i(q, tau_prime, prec);
    acb_modular_theta_sum(t1, t2, t3, t4, w, 1, q, 1, prec);

    /* (t2 t3 t4) ^ 8 * q^2 */
    acb_mul(t1, t2, t3, prec);
    acb_mul(t1, t1, t4, prec);
    acb_mul(t1, t1, t1, prec);
    acb_mul(t1, t1, t1, prec);
    acb_mul(t1, t1, q, prec);
    acb_mul(t1, t1, t1, prec);
    acb_mul_2exp_si(t1, t1, -8);

    if (!fmpz_is_zero(&g->c))
    {
        acb_mul_fmpz(t2, tau, &g->c, prec);
        acb_add_fmpz(t2, t2, &g->d, prec);
        acb_pow_ui(t2, t2, 12, prec);
        acb_div(t1, t1, t2, prec);
    }

    acb_set(z, t1);

    psl2z_clear(g);
    arf_clear(one_minus_eps);
    acb_clear(tau_prime);
    acb_clear(t1);
    acb_clear(t2);
    acb_clear(t3);
    acb_clear(t4);
    acb_clear(w);
    acb_clear(q);
}

