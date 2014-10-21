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
#include "acb_poly.h"

void
acb_modular_theta_zpx_notransform(acb_ptr theta1, acb_ptr theta2,
    acb_ptr theta3, acb_ptr theta4, const acb_t z, const acb_t tau,
    long len, long prec)
{
    acb_t q, q4, w;
    int w_is_unit;

    acb_init(q);
    acb_init(q4);
    acb_init(w);

    /* compute q_{1/4}, q */
    acb_mul_2exp_si(q4, tau, -2);
    acb_exp_pi_i(q4, q4, prec);
    acb_pow_ui(q, q4, 4, prec);

    /* compute w */
    acb_exp_pi_i(w, z, prec);
    w_is_unit = arb_is_zero(acb_imagref(z));

    /* evaluate theta functions */
    acb_modular_theta_sum(theta1, theta2, theta3, theta4,
        w, w_is_unit, q, len, prec);
    _acb_vec_scalar_mul(theta1, theta1, len, q4, prec);
    _acb_vec_scalar_mul(theta2, theta2, len, q4, prec);

    acb_clear(q);
    acb_clear(q4);
    acb_clear(w);
}

void
acb_modular_elliptic_p_zpx(acb_ptr r, const acb_t z, const acb_t tau, long len, long prec)
{
    acb_t t01, t02, t03, t04;
    acb_ptr tz1, tz2, tz3, tz4;
    acb_t t;

    if (len < 1)
        return;

    if (len == 1)
    {
        acb_modular_elliptic_p(r, z, tau, prec);
        return;
    }

    acb_init(t);

    acb_init(t01);
    acb_init(t02);
    acb_init(t03);
    acb_init(t04);

    tz1 = _acb_vec_init(len);
    tz2 = _acb_vec_init(len);
    tz3 = _acb_vec_init(len);
    tz4 = _acb_vec_init(len);

    acb_modular_theta_zpx_notransform(tz1, tz2, tz3, tz4, z, tau, len, prec);

    /* [theta_4(z) / theta_1(z)]^2 */
    _acb_poly_div_series(tz2, tz4, len, tz1, len, len, prec);
    _acb_poly_mullow(tz1, tz2, len, tz2, len, len, prec);

    acb_zero(t);
    acb_modular_theta_notransform(t01, t02, t03, t04, t, tau, prec);

    /* [theta_2(0) * theta_3(0)] ^2 */
    acb_mul(t, t02, t03, prec);
    acb_mul(t, t, t, prec);
    _acb_vec_scalar_mul(tz1, tz1, len, t, prec);

    /* - [theta_2(0)^4 + theta_3(0)^4] / 3 */
    acb_pow_ui(t02, t02, 4, prec);
    acb_pow_ui(t03, t03, 4, prec);
    acb_add(t, t02, t03, prec);
    acb_div_ui(t, t, 3, prec);
    acb_sub(tz1, tz1, t, prec);

    /* times pi^2 */
    acb_const_pi(t, prec);
    acb_mul(t, t, t, prec);
    _acb_vec_scalar_mul(r, tz1, len, t, prec);

    acb_clear(t);

    acb_clear(t01);
    acb_clear(t02);
    acb_clear(t03);
    acb_clear(t04);

    _acb_vec_clear(tz1, len);
    _acb_vec_clear(tz2, len);
    _acb_vec_clear(tz3, len);
    _acb_vec_clear(tz4, len);
}

