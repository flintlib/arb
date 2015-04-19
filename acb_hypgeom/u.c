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
acb_hypgeom_u_1f1_series(acb_poly_t res,
    const acb_poly_t a, const acb_poly_t b, const acb_poly_t z,
    long len, long prec)
{
    acb_poly_t s, u, A, B;
    acb_poly_struct aa[3];
    arb_t c;
    long wlen;
    int singular;

    acb_poly_init(s);
    acb_poly_init(u);
    acb_poly_init(A);
    acb_poly_init(B);
    acb_poly_init(aa + 0);
    acb_poly_init(aa + 1);
    acb_poly_init(aa + 2);
    arb_init(c);

    singular = (b->length == 0) || acb_is_int(b->coeffs);
    wlen = len + (singular != 0);

    /* A = rgamma(a-b+1) * 1F~1(a,b,z) */
    acb_poly_sub(u, a, b, prec);
    acb_poly_add_si(u, u, 1, prec);
    acb_poly_rgamma_series(A, u, wlen, prec);

    /* todo: handle a = 1 efficiently */
    acb_poly_set(aa, a);
    acb_poly_set(aa + 1, b);
    acb_poly_one(aa + 2);
    acb_hypgeom_pfq_series_direct(s, aa, 1, aa + 1, 2, z, 1, -1, wlen, prec);
    acb_poly_mullow(A, A, s, wlen, prec);

    /* B = rgamma(a) * 1F~1(a-b+1,2-b,z) * z^(1-b) */
    acb_poly_set(aa, u);
    acb_poly_add_si(aa + 1, b, -2, prec);
    acb_poly_neg(aa + 1, aa + 1);
    acb_hypgeom_pfq_series_direct(s, aa, 1, aa + 1, 2, z, 1, -1, wlen, prec);
    acb_poly_rgamma_series(B, a, wlen, prec);
    acb_poly_mullow(B, B, s, wlen, prec);

    acb_poly_add_si(u, b, -1, prec);
    acb_poly_neg(u, u);
    acb_poly_pow_series(s, z, u, wlen, prec);
    acb_poly_mullow(B, B, s, wlen, prec);

    acb_poly_sub(A, A, B, prec);

    /* multiply by pi csc(pi b) */
    acb_poly_sin_pi_series(B, b, wlen, prec);

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
    acb_poly_clear(aa + 0);
    acb_poly_clear(aa + 1);
    acb_poly_clear(aa + 2);
    arb_clear(c);
}

void
acb_hypgeom_u_1f1(acb_t res, const acb_t a, const acb_t b, const acb_t z, long prec)
{
    if (acb_is_int(b))
    {
        acb_poly_t aa, bb, zz;

        acb_poly_init(aa);
        acb_poly_init(bb);
        acb_poly_init(zz);

        acb_poly_set_acb(aa, a);
        acb_poly_set_coeff_acb(bb, 0, b);
        acb_poly_set_coeff_si(bb, 1, 1);
        acb_poly_set_acb(zz, z);

        acb_hypgeom_u_1f1_series(zz, aa, bb, zz, 1, prec);
        acb_poly_get_coeff_acb(res, zz, 0);

        acb_poly_clear(aa);
        acb_poly_clear(bb);
        acb_poly_clear(zz);
    }
    else
    {
        acb_t t, u, v;
        acb_struct aa[3];

        acb_init(t);
        acb_init(u);
        acb_init(v);
        acb_init(aa + 0);
        acb_init(aa + 1);
        acb_init(aa + 2);

        acb_set(aa, a);
        acb_set(aa + 1, b);
        acb_one(aa + 2);
        acb_hypgeom_pfq_direct(u, aa, 1, aa + 1, 2, z, -1, prec);
        acb_sub(aa, a, b, prec);
        acb_add_ui(aa, aa, 1, prec);
        acb_sub_ui(aa + 1, b, 2, prec);
        acb_neg(aa + 1, aa + 1);
        acb_hypgeom_pfq_direct(v, aa, 1, aa + 1, 2, z, -1, prec);

        acb_sub_ui(aa + 1, b, 1, prec);

        /* rgamma(a-b+1) * gamma(1-b) * u */
        acb_rgamma(t, aa, prec);
        acb_mul(u, u, t, prec);
        acb_neg(t, aa + 1);
        acb_gamma(t, t, prec);
        acb_mul(u, u, t, prec);

        /* rgamma(a) * gamma(b-1) * z^(1-b) * v */
        acb_rgamma(t, a, prec);
        acb_mul(v, v, t, prec);
        acb_gamma(t, aa + 1, prec);
        acb_mul(v, v, t, prec);
        acb_neg(t, aa + 1);
        acb_pow(t, z, t, prec);
        acb_mul(v, v, t, prec);

        acb_add(res, u, v, prec);

        acb_clear(t);
        acb_clear(u);
        acb_clear(v);
        acb_clear(aa + 0);
        acb_clear(aa + 1);
        acb_clear(aa + 2);
    }
}

void
acb_hypgeom_u(acb_t res, const acb_t a, const acb_t b, const acb_t z, long prec)
{
    acb_t t;
    acb_init(t);

    acb_sub(t, a, b, prec);
    acb_add_ui(t, t, 1, prec);

    if ((acb_is_int(a) && arf_sgn(arb_midref(acb_realref(a))) <= 0) ||
        (acb_is_int(t) && arf_sgn(arb_midref(acb_realref(t))) <= 0) ||
        acb_hypgeom_u_use_asymp(z, prec))
    {
        acb_neg(t, a);
        acb_pow(t, z, t, prec);
        acb_hypgeom_u_asymp(res, a, b, z, -1, prec);
        acb_mul(res, res, t, prec);
    }
    else
    {
        acb_hypgeom_u_1f1(res, a, b, z, prec);
    }

    acb_clear(t);
}

