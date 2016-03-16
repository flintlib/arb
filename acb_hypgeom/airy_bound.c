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
arb_bound_exp_neg(mag_t b, const arb_t x)
{
    arb_t t;
    arb_init(t);
    arf_set_mag(arb_midref(t), arb_radref(x));
    arf_sub(arb_midref(t), arb_midref(x), arb_midref(t), MAG_BITS, ARF_RND_FLOOR);
    arf_neg(arb_midref(t), arb_midref(t));
    arb_exp(t, t, MAG_BITS);
    arb_get_mag(b, t);
    arb_clear(t);
}

/* todo -- should be lt in asymp code? */
static int
arg_le_2pi3(const acb_t z, const acb_t zeta)
{
    if (arb_is_nonnegative(acb_realref(z)))
        return 1;

    if (arb_is_positive(acb_imagref(z)) &&
        arb_is_nonnegative(acb_imagref(zeta)))
        return 1;

    if (arb_is_negative(acb_imagref(z)) &&
        arb_is_nonpositive(acb_imagref(zeta)))
        return 1;

    return 0;
}

void
acb_hypgeom_airy_bound_9_7_17(mag_t bound, const acb_t z, const acb_t zeta)
{
    mag_t D, t, u, v, zeta_lower;

    mag_init(D);
    mag_init(t);
    mag_init(u);
    mag_init(v);
    mag_init(zeta_lower);

    acb_get_mag_lower(zeta_lower, zeta);

    /* 2 chi(1) exp(7 pi / (72 |zeta|)) * c_1 */
    /* simplified bound */
    if (mag_cmp_2exp_si(zeta_lower, -1) >= 0)
        mag_one(D);
    else
        mag_inf(D);

    if (!arg_le_2pi3(z, zeta))
    {
        arb_get_mag_lower(u, acb_realref(zeta));
        arb_get_mag(v, acb_imagref(zeta));

        /* exp(7 pi / (36 u)) < exp(5/(8 u)) */
        mag_set_ui_2exp_si(t, 5, -3);
        mag_div(t, t, u);
        mag_exp(t, t);

        /* |1/cos(arg(zeta))| = sqrt(1+(v/u)^2) */
        mag_div(v, v, u);
        mag_mul(v, v, v);
        mag_one(u);
        mag_add(v, v, u);
        mag_sqrt(v, v);
        mag_mul(t, t, v);
        /* c_1 * 4 chi(1) < 0.62 < 1 -- do nothing */

        mag_max(D, D, t);
    }

    /* exp(-zeta) / (2 sqrt(pi)) * (1 + D / |zeta|) */

    /* t = exp(-zeta) / (2 sqrt(pi)) < exp(-zeta) * (73/256) */
    arb_bound_exp_neg(t, acb_realref(zeta));
    mag_mul_ui(t, t, 73);
    mag_mul_2exp_si(t, t, -8);

    /* u = 1 + D / |zeta| */
    mag_div(u, D, zeta_lower);
    mag_one(v);
    mag_add(u, u, v);

    mag_mul(bound, t, u);

    mag_clear(D);
    mag_clear(t);
    mag_clear(u);
    mag_clear(v);
    mag_clear(zeta_lower);
}

void
acb_hypgeom_airy_bound(mag_t ai, mag_t aip, mag_t bi, mag_t bip, const acb_t z)
{
    acb_t zeta;
    acb_t z1, z2;
    arf_srcptr zre, zim;
    mag_t A, B, D, zlo, zhi;
    slong wp;

    if (acb_contains_zero(z))
    {
        if (ai != NULL) mag_inf(ai);
        if (aip != NULL) mag_inf(aip);
        if (bi != NULL) mag_inf(bi);
        if (bip != NULL) mag_inf(bip);
        return;
    }

    acb_init(zeta);
    mag_init(A);
    mag_init(B);
    mag_init(D);
    mag_init(zlo);
    mag_init(zhi);

    wp = MAG_BITS * 2;
    zre = arb_midref(acb_realref(z));
    zim = arb_midref(acb_imagref(z));

    /* near the negative half line, use 
       |Ai(-z)|, |Bi(-z)| <= |Ai(z exp(pi i/3))| + |Ai(z exp(-pi i/3))| */
    if (arf_sgn(zre) < 0 && arf_cmpabs(zre, zim) > 0)
    {
        acb_init(z1);
        acb_init(z2);

        arb_sqrt_ui(acb_imagref(z1), 3, wp);
        arb_one(acb_realref(z1));
        acb_mul_2exp_si(z1, z1, -1);
        acb_conj(z2, z1);

        acb_neg_round(zeta, z, wp);
        acb_mul(z1, z1, zeta, wp);

        acb_sqrt(zeta, zeta, wp);
        acb_cube(zeta, zeta, wp);
        acb_mul_2exp_si(zeta, zeta, 1);
        acb_div_ui(zeta, zeta, 3, wp);
        acb_mul_onei(zeta, zeta);

        acb_hypgeom_airy_bound_9_7_17(A, z1, zeta);

        /* conjugate symmetry */
        if (arb_is_zero(acb_imagref(z)))
        {
            mag_mul_2exp_si(A, A, 1);
        }
        else
        {
            acb_mul(z2, z2, zeta, wp);
            acb_neg(zeta, zeta);
            acb_hypgeom_airy_bound_9_7_17(D, z2, zeta);
            mag_add(A, A, D);
        }

        mag_set(B, A);

        acb_clear(z1);
        acb_clear(z2);
    }
    else
    {
        acb_set_round(zeta, z, wp);
        acb_sqrt(zeta, zeta, wp);
        acb_cube(zeta, zeta, wp);
        acb_mul_2exp_si(zeta, zeta, 1);
        acb_div_ui(zeta, zeta, 3, wp);

        acb_hypgeom_airy_bound_9_7_17(A, z, zeta);

        /* Use Bi(z) = w_1 Ai(z) + 2 w_2 Ai(z exp(+/- 2pi i / 3)),
           where w_1, w_2 are roots of unity */
        if (bi != NULL || bip != NULL)
        {
            acb_init(z1);

            arb_sqrt_ui(acb_imagref(z1), 3, wp);
            arb_set_si(acb_realref(z1), -1);
            acb_mul_2exp_si(z1, z1, -1);

            /* multiply by exp(-2 pi i / 3) in upper half plane
               and by exp(2 pi i / 3) in lower half plane, to stay close
               to positive reals */
            if (arf_sgn(zim) >= 0)
                acb_conj(z1, z1);

            acb_mul(z1, z1, z, wp);
            acb_neg(zeta, zeta);   /* same effect regardless of exp(+/-2 pi i/3) */

            acb_hypgeom_airy_bound_9_7_17(B, z1, zeta);
            mag_mul_2exp_si(B, B, 1);
            mag_add(B, B, A);

            acb_clear(z1);
        }
    }

    acb_get_mag(zhi, z);
    acb_get_mag_lower(zlo, z);

    /* bound |z|^(1/4) */
    mag_sqrt(zhi, zhi);
    mag_sqrt(zhi, zhi);

    /* bound |z|^(-1/4) */
    mag_rsqrt(zlo, zlo);
    mag_sqrt(zlo, zlo);

    if (ai != NULL) mag_mul(ai, A, zlo);
    if (aip != NULL) mag_mul(aip, A, zhi);
    if (bi != NULL) mag_mul(bi, B, zlo);
    if (bip != NULL) mag_mul(bip, B, zhi);

    acb_clear(zeta);
    mag_clear(A);
    mag_clear(B);
    mag_clear(D);
    mag_clear(zlo);
    mag_clear(zhi);
}

