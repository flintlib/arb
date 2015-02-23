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

void
acb_hypgeom_bessel_k_asymp(acb_t res, const acb_t nu, const acb_t z, long prec)
{
    acb_t t, a, b, w;

    acb_init(t);
    acb_init(a);
    acb_init(b);
    acb_init(w);

    acb_one(a);
    acb_mul_2exp_si(a, a, -1);
    acb_add(a, a, nu, prec);

    acb_mul_2exp_si(b, nu, 1);
    acb_add_ui(b, b, 1, prec);

    acb_mul_2exp_si(w, z, 1);

    acb_hypgeom_u_asymp(t, a, b, w, -1, prec);

    acb_neg(w, z);
    acb_exp(w, w, prec);    
    acb_mul(t, t, w, prec);

    acb_const_pi(w, prec);
    acb_div(w, w, z, prec);
    acb_mul_2exp_si(w, w, -1);
    acb_sqrt(w, w, prec);

    acb_mul(res, t, w, prec);

    acb_clear(t);
    acb_clear(a);
    acb_clear(b);
    acb_clear(w);
}

void
acb_hypgeom_bessel_k_0f1_series(acb_poly_t res,
    const acb_poly_t nu, const acb_poly_t z,
    long len, long prec)
{
    acb_poly_t s, u, A, B;
    acb_poly_struct b[2];
    arb_t c;
    long wlen;
    int singular;

    acb_poly_init(s);
    acb_poly_init(u);
    acb_poly_init(A);
    acb_poly_init(B);
    acb_poly_init(b + 0);
    acb_poly_init(b + 1);
    arb_init(c);

    singular = (nu->length == 0) || acb_is_int(nu->coeffs);
    wlen = len + (singular != 0);

    /* A = (z/2)^nu, B = 1/A */
    acb_poly_scalar_mul_2exp_si(A, z, -1);
    acb_poly_pow_series(A, A, nu, wlen, prec);
    acb_poly_inv_series(B, A, wlen, prec);

    acb_poly_mullow(u, z, z, wlen, prec);
    acb_poly_scalar_mul_2exp_si(u, u, -2);

    acb_poly_one(b + 1);
    acb_poly_add_si(b + 0, nu, 1, prec);
    acb_hypgeom_pfq_series_direct(s, NULL, 0, b, 2, u, -1, wlen, prec, 1);
    acb_poly_mullow(A, A, s, wlen, prec);

    acb_poly_add_si(b + 0, nu, -1, prec);
    acb_poly_neg(b + 0, b + 0);
    acb_hypgeom_pfq_series_direct(s, NULL, 0, b, 2, u, -1, wlen, prec, 1);
    acb_poly_mullow(B, B, s, wlen, prec);

    acb_poly_sub(A, B, A, prec);
    acb_poly_scalar_mul_2exp_si(A, A, -1);

    /* multiply by pi csc(pi nu) */
    acb_poly_sin_pi_series(B, nu, wlen, prec);

    if (singular)
    {
        acb_poly_shift_right(A, A, 1);
        acb_poly_shift_right(B, B, 1);
    }

    acb_poly_div_series(res, A, B, len, prec);

    arb_const_pi(c, prec);
    _acb_vec_scalar_mul_arb(res->coeffs, res->coeffs, res->length, c, prec);

    acb_poly_clear(s);
    acb_poly_clear(u);
    acb_poly_clear(A);
    acb_poly_clear(B);
    acb_poly_clear(b + 0);
    acb_poly_clear(b + 1);
    arb_clear(c);
}

void
acb_hypgeom_bessel_k_0f1(acb_t res, const acb_t nu, const acb_t z, long prec)
{
    if (acb_is_int(nu))
    {
        acb_poly_t nux, zx, rx;

        acb_poly_init(nux);
        acb_poly_init(zx);
        acb_poly_init(rx);

        acb_poly_set_coeff_acb(nux, 0, nu);
        acb_poly_set_coeff_si(nux, 1, 1);
        acb_poly_set_acb(zx, z);

        acb_hypgeom_bessel_k_0f1_series(rx, nux, zx, 1, prec);

        acb_poly_get_coeff_acb(res, rx, 0);

        acb_poly_clear(nux);
        acb_poly_clear(zx);
        acb_poly_clear(rx);
    }
    else
    {
        acb_t t, u, v, w;
        acb_struct b[2];

        acb_init(t);
        acb_init(u);
        acb_init(v);
        acb_init(w);
        acb_init(b + 0);
        acb_init(b + 1);

        /* u = 0F1(1+nu), v = 0F1(1-nu) */
        acb_mul(t, z, z, prec);
        acb_mul_2exp_si(t, t, -2);
        acb_add_ui(b, nu, 1, prec);
        acb_one(b + 1);
        acb_hypgeom_pfq_direct(u, NULL, 0, b, 2, t, -1, prec);
        acb_sub_ui(b, nu, 1, prec);
        acb_neg(b, b);
        acb_hypgeom_pfq_direct(v, NULL, 0, b, 2, t, -1, prec);

        /* v = v * gamma(nu) / (z/2)^nu */
        acb_mul_2exp_si(t, z, -1);
        acb_pow(t, t, nu, prec);
        acb_gamma(w, nu, prec);
        acb_mul(v, v, w, prec);
        acb_div(v, v, t, prec);

        /* u = u * t * pi / (gamma(nu) * nu * sin(pi nu)) */
        acb_mul(u, u, t, prec);
        acb_const_pi(t, prec);
        acb_mul(u, u, t, prec);
        acb_sin_pi(t, nu, prec);
        acb_mul(t, t, w, prec);
        acb_mul(t, t, nu, prec);
        acb_div(u, u, t, prec);

        acb_sub(res, v, u, prec);
        acb_mul_2exp_si(res, res, -1);

        acb_clear(t);
        acb_clear(u);
        acb_clear(v);
        acb_clear(w);
        acb_clear(b + 0);
        acb_clear(b + 1);
    }
}

void
acb_hypgeom_bessel_k(acb_t res, const acb_t nu, const acb_t z, long prec)
{
    mag_t zmag;

    mag_init(zmag);
    acb_get_mag(zmag, z);

    if (mag_cmp_2exp_si(zmag, 4) < 0 ||
        (mag_cmp_2exp_si(zmag, 64) < 0 && 2 * mag_get_d(zmag) < prec))
        acb_hypgeom_bessel_k_0f1(res, nu, z, prec);
    else
        acb_hypgeom_bessel_k_asymp(res, nu, z, prec);

    mag_clear(zmag);
}

