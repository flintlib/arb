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

#include "acb_hypgeom.h"

/* todo: negative z and integer nu */
/* todo: integer nu? */
/* todo: integer nu+1/2 */
/* todo: change precision when one term is small */

void
acb_hypgeom_bessel_j_asymp(acb_t res, const acb_t nu, const acb_t z, long prec)
{
    acb_t A1, A2, B1, B2, s, t, u;
    int is_real;

    acb_init(A1);
    acb_init(A2);
    acb_init(B1);
    acb_init(B2);
    acb_init(s);
    acb_init(t);
    acb_init(u);

    is_real = acb_is_real(nu) && acb_is_real(z)
        && (acb_is_int(nu) || arb_is_positive(acb_realref(z)));

    /* s = 1/2 + nu */
    acb_one(s);
    acb_mul_2exp_si(s, s, -1);
    acb_add(s, s, nu, prec);

    /* t = 1 + 2 nu */
    acb_mul_2exp_si(t, nu, 1);
    acb_add_ui(t, t, 1, prec);

    acb_mul_onei(u, z);
    acb_mul_2exp_si(u, u, 1);
    acb_hypgeom_u_asymp(B2, s, t, u, -1, prec);
    acb_neg(u, u);
    acb_hypgeom_u_asymp(B1, s, t, u, -1, prec);

    if (arb_is_positive(acb_realref(z)))
    {
        /* -(2nu+1)/4 pi + z */
        acb_mul_2exp_si(t, nu, 1);
        acb_add_ui(t, t, 1, prec);
        acb_mul_2exp_si(t, t, -2);
        acb_neg(t, t);
        acb_const_pi(u, prec);
        acb_mul(t, t, u, prec);
        acb_add(t, t, z, prec);
        acb_mul_onei(t, t);

        acb_exp(A1, t, prec);

        if (acb_is_real(nu) && acb_is_real(z))
        {
            acb_conj(A2, A1);
        }
        else
        {
            acb_neg(t, t);
            acb_exp(A2, t, prec);
        }

        acb_mul(A1, A1, B1, prec);
        acb_mul(A2, A2, B2, prec);
        acb_add(A1, A1, A2, prec);

        /* divide by sqrt(2 pi z) */
        acb_const_pi(u, prec);
        acb_mul_2exp_si(u, u, 1);
        acb_mul(u, u, z, prec);
        acb_rsqrt(u, u, prec);
        acb_mul(A1, A1, u, prec);

        acb_set(res, A1);
    }
    else
    {
        /* general case */

        /* exp(iz), exp(-iz) */
        acb_mul_onei(t, z);
        acb_exp(u, t, prec);
        acb_mul(B1, B1, u, prec);

        if (acb_is_real(z))
        {
            acb_conj(u, u);
        }
        else
        {
            acb_neg(t, t);
            acb_exp(u, t, prec);
        }

        acb_mul(B2, B2, u, prec);

        /* s = -(1/2+nu) */
        acb_one(s);
        acb_mul_2exp_si(s, s, -1);
        acb_add(s, s, nu, prec);
        acb_neg(s, s);
        acb_mul_onei(t, z);
        acb_pow(A1, t, s, prec);
        acb_neg(t, t);
        acb_pow(A2, t, s, prec);

        acb_mul(A1, A1, B1, prec);
        acb_mul(A2, A2, B2, prec);
        acb_add(A1, A1, A2, prec);

        acb_pow(t, z, nu, prec);
        acb_mul(A1, A1, t, prec);

        /* divide by sqrt(2pi) */
        acb_const_pi(u, prec);
        acb_mul_2exp_si(u, u, 1);
        acb_rsqrt(u, u, prec);
        acb_mul(A1, A1, u, prec);

        acb_set(res, A1);
    }

    if (is_real)
        arb_zero(acb_imagref(res));

    acb_clear(A1);
    acb_clear(A2);
    acb_clear(B1);
    acb_clear(B2);
    acb_clear(s);
    acb_clear(t);
    acb_clear(u);
}

void
acb_hypgeom_bessel_j_0f1(acb_t res, const acb_t nu, const acb_t z, long prec)
{
    acb_struct b[2];
    acb_t w, c, t;

    if (acb_is_int(nu) && arb_is_negative(acb_realref(nu)))
    {
        acb_init(t);
        acb_neg(t, nu);

        acb_hypgeom_bessel_j_0f1(res, t, z, prec);

        acb_mul_2exp_si(t, t, -1);
        if (!acb_is_int(t))
            acb_neg(res, res);

        acb_clear(t);
        return;
    }

    acb_init(b + 0);
    acb_init(b + 1);
    acb_init(w);
    acb_init(c);
    acb_init(t);

    acb_add_ui(b + 0, nu, 1, prec);
    acb_one(b + 1);

    /* (z/2)^nu / gamma(nu+1) */
    acb_mul_2exp_si(c, z, -1);
    acb_pow(c, c, nu, prec);
    acb_rgamma(t, b + 0, prec);
    acb_mul(c, t, c, prec);

    /* -z^2/4 */
    acb_mul(w, z, z, prec);
    acb_mul_2exp_si(w, w, -2);
    acb_neg(w, w);

    acb_hypgeom_pfq_direct(t, NULL, 0, b, 2, w, -1, prec);

    acb_mul(res, t, c, prec);

    acb_clear(b + 0);
    acb_clear(b + 1);
    acb_clear(w);
    acb_clear(c);
    acb_clear(t);
}

/*

The asymptotic series can be used roughly when

[(1+log(2))/log(2) = 2.44269504088896] * z > p

We are a bit more conservative and use the factor 2.
*/

void
acb_hypgeom_bessel_j(acb_t res, const acb_t nu, const acb_t z, long prec)
{
    mag_t zmag;

    mag_init(zmag);
    acb_get_mag(zmag, z);

    if (mag_cmp_2exp_si(zmag, 4) < 0 ||
        (mag_cmp_2exp_si(zmag, 64) < 0 && 2 * mag_get_d(zmag) < prec))
        acb_hypgeom_bessel_j_0f1(res, nu, z, prec);
    else
        acb_hypgeom_bessel_j_asymp(res, nu, z, prec);

    mag_clear(zmag);
}

