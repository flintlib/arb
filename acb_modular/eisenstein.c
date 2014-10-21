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
acb_modular_eisenstein(acb_ptr r, const acb_t tau, long len, long prec)
{
    psl2z_t g;
    arf_t one_minus_eps;
    acb_t tau_prime, t1, t2, t3, t4, w, q;
    long m, n;

    if (len < 1)
        return;

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

    /* fourth powers of the theta functions (a, b, c) */
    acb_mul(t2, t2, t2, prec);
    acb_mul(t2, t2, t2, prec);
    acb_mul(t2, t2, q, prec);

    acb_mul(t3, t3, t3, prec);
    acb_mul(t3, t3, t3, prec);

    acb_mul(t4, t4, t4, prec);
    acb_mul(t4, t4, t4, prec);

    /* c2 = pi^4 * (a^8 + b^8 + c^8) / 30 */
    /* c3 = pi^6 * (b^12 + c^12 - 3a^8 * (b^4+c^4)) / 180 */

    /* r = a^8 */
    acb_mul(r, t2, t2, prec);

    if (len > 1)
    {
        /* r[1] = -3 a^8 * (b^4 + c^4) */
        acb_add(r + 1, t3, t4, prec);
        acb_mul(r + 1, r + 1, r, prec);
        acb_mul_si(r + 1, r + 1, -3, prec);
    }

    /* b^8 */
    acb_mul(t1, t3, t3, prec);
    acb_add(r, r, t1, prec);

    /* b^12 */
    if (len > 1)
        acb_addmul(r + 1, t1, t3, prec);

    /* c^8 */
    acb_mul(t1, t4, t4, prec);
    acb_add(r, r, t1, prec);

    /* c^12 */
    if (len > 1)
        acb_addmul(r + 1, t1, t4, prec);

    acb_const_pi(t1, prec);
    acb_mul(t1, t1, t1, prec);
    acb_mul(t2, t1, t1, prec);
    acb_mul(r, r, t2, prec);
    acb_div_ui(r, r, 30, prec);

    if (len > 1)
    {
        acb_mul(t2, t2, t1, prec);
        acb_mul(r + 1, r + 1, t2, prec);
        acb_div_ui(r + 1, r + 1, 189, prec);
    }

    /* apply modular transformation */
    if (!fmpz_is_zero(&g->c))
    {

        acb_mul_fmpz(t1, tau, &g->c, prec);
        acb_add_fmpz(t1, t1, &g->d, prec);
        acb_inv(t1, t1, prec);
        acb_mul(t1, t1, t1, prec);
        acb_mul(t2, t1, t1, prec);
        acb_mul(r, r, t2, prec);

        if (len > 1)
        {
            acb_mul(t2, t1, t2, prec);
            acb_mul(r + 1, r + 1, t2, prec);
        }
    }

    /* compute more coefficients using recurrence */
    for (n = 4; n < len + 2; n++)
    {
        acb_zero(r + n - 2);

        m = 2;
        for (m = 2; m * 2 < n; m++)
            acb_addmul(r + n - 2, r + m - 2, r + n - m - 2, prec);

        acb_mul_2exp_si(r + n - 2, r + n - 2, 1);

        if (n % 2 == 0)
            acb_addmul(r + n - 2, r + n / 2 - 2, r + n / 2 - 2, prec);

        acb_mul_ui(r + n - 2, r + n - 2, 3, prec);
        acb_div_ui(r + n - 2, r + n - 2, (2 * n + 1) * (n - 3), prec);
    }

    /* convert c's to G's */
    for (n = 0; n < len; n++)
        acb_div_ui(r + n, r + n, 2 * n + 3, prec);

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

