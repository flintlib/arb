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
_acb_hypgeom_2f1_inf_limit(acb_t res, const acb_poly_t a, const acb_poly_t b,
    const acb_poly_t c, const acb_poly_t z, long prec)
{
    acb_poly_t ba, ca, cb, ac1, bc1, ab1, ba1, w, t, u, v;

    acb_poly_init(ba);
    acb_poly_init(ca); acb_poly_init(cb);
    acb_poly_init(ac1); acb_poly_init(bc1);
    acb_poly_init(ab1); acb_poly_init(ba1);
    acb_poly_init(w); acb_poly_init(t);
    acb_poly_init(u); acb_poly_init(v);

    acb_poly_sub(ba, b, a, prec); /* b - a */
    acb_poly_sub(ca, c, a, prec); /* c - a */
    acb_poly_sub(cb, c, b, prec); /* c - b */

    acb_poly_add_si(ac1, ca, -1, prec); acb_poly_neg(ac1, ac1); /* a - c + 1 */
    acb_poly_add_si(bc1, cb, -1, prec); acb_poly_neg(bc1, bc1); /* b - c + 1 */
    acb_poly_add_si(ab1, ba, -1, prec); acb_poly_neg(ab1, ab1); /* a - b + 1 */
    acb_poly_add_si(ba1, ba, 1, prec);                          /* b - a + 1 */

    acb_poly_inv_series(w, z, 2, prec);

    acb_hypgeom_2f1_series_direct(t, a, ac1, ab1, w, 1, 2, prec);

    acb_poly_rgamma_series(v, b, 2, prec);
    acb_poly_mullow(t, t, v, 2, prec);
    acb_poly_rgamma_series(v, ca, 2, prec);
    acb_poly_mullow(t, t, v, 2, prec);
    acb_poly_neg(v, z);
    acb_poly_neg(ca, a);
    acb_poly_pow_series(v, v, ca, 2, prec);
    acb_poly_mullow(t, t, v, 2, prec);

    acb_hypgeom_2f1_series_direct(u, b, bc1, ba1, w, 1, 2, prec);

    acb_poly_rgamma_series(v, a, 2, prec);
    acb_poly_mullow(u, u, v, 2, prec);
    acb_poly_rgamma_series(v, cb, 2, prec);
    acb_poly_mullow(u, u, v, 2, prec);
    acb_poly_neg(v, z);
    acb_poly_neg(cb, b);
    acb_poly_pow_series(v, v, cb, 2, prec);
    acb_poly_mullow(u, u, v, 2, prec);

    acb_poly_sub(t, t, u, prec);

    /* todo: this is just multiplying by +/- 1 when
       not computing any further derivatives */
    acb_poly_sin_pi_series(v, ba, 2, prec);

    {
        acb_t tt;
        acb_init(tt);
        acb_poly_get_coeff_acb(tt, t, 1);
        acb_poly_get_coeff_acb(res, v, 1);
        acb_div(res, tt, res, prec);
        acb_const_pi(tt, prec);
        acb_mul(res, res, tt, prec);
        acb_clear(tt);
    }

    acb_poly_clear(ba);
    acb_poly_clear(ca); acb_poly_clear(cb);
    acb_poly_clear(ac1); acb_poly_clear(bc1);
    acb_poly_clear(ab1); acb_poly_clear(ba1);
    acb_poly_clear(w); acb_poly_clear(t);
    acb_poly_clear(u); acb_poly_clear(v);
}

void
acb_hypgeom_2f1_inf_limit(acb_t res, const acb_t a, const acb_t b,
    const acb_t c, const acb_t z, int regularized, long prec)
{
    acb_poly_t aa, bb, cc, zz;

    if (acb_contains_zero(z))
    {
        acb_indeterminate(res);
        return;
    }

    if (!regularized)
    {
        acb_t t;
        acb_init(t);
        acb_gamma(t, c, prec);
        acb_hypgeom_2f1_inf_limit(res, a, b, c, z, 1, prec);
        acb_mul(res, res, t, prec);
        acb_clear(t);
        return;
    }

    acb_poly_init(aa);
    acb_poly_init(bb);
    acb_poly_init(cc);
    acb_poly_init(zz);

    acb_poly_set_acb(aa, a);
    acb_poly_set_acb(bb, b);
    acb_poly_set_acb(cc, c);
    acb_poly_set_acb(zz, z);

    acb_poly_set_coeff_si(aa, 1, 1);

    _acb_hypgeom_2f1_inf_limit(res, aa, bb, cc, zz, prec);

    acb_poly_clear(aa);
    acb_poly_clear(bb);
    acb_poly_clear(cc);
    acb_poly_clear(zz);
}

void
acb_hypgeom_2f1_inf(acb_t res, const acb_t a, const acb_t b,
    const acb_t c, const acb_t z, int regularized, long prec)
{
    acb_t ba, ca, cb, ac1, bc1, ab1, ba1, w, t, u, v;

    if (acb_contains_zero(z))
    {
        acb_indeterminate(res);
        return;
    }

    if (!regularized)
    {
        acb_init(t);
        acb_gamma(t, c, prec);
        acb_hypgeom_2f1_inf(res, a, b, c, z, 1, prec);
        acb_mul(res, res, t, prec);
        acb_clear(t);
        return;
    }

    acb_init(ba);
    acb_sub(ba, b, a, prec); /* b - a */

    if (acb_is_int(ba))
    {
        acb_hypgeom_2f1_inf_limit(res, a, b, c, z, regularized, prec);
    }
    else
    {
        acb_init(ca); acb_init(cb);
        acb_init(ac1); acb_init(bc1);
        acb_init(ab1); acb_init(ba1);
        acb_init(w); acb_init(t);
        acb_init(u); acb_init(v);

        acb_sub(ca, c, a, prec); /* c - a */
        acb_sub(cb, c, b, prec); /* c - b */

        acb_sub_ui(ac1, ca, 1, prec); acb_neg(ac1, ac1); /* a - c + 1 */
        acb_sub_ui(bc1, cb, 1, prec); acb_neg(bc1, bc1); /* b - c + 1 */
        acb_sub_ui(ab1, ba, 1, prec); acb_neg(ab1, ab1); /* a - b + 1 */
        acb_add_ui(ba1, ba, 1, prec);                    /* b - a + 1 */

        acb_inv(w, z, prec);

        acb_hypgeom_2f1_direct(t, a, ac1, ab1, w, 1, prec);

        acb_rgamma(v, b, prec);
        acb_mul(t, t, v, prec);
        acb_rgamma(v, ca, prec);
        acb_mul(t, t, v, prec);
        acb_neg(v, z);
        acb_neg(ca, a);
        acb_pow(v, v, ca, prec);
        acb_mul(t, t, v, prec);

        acb_hypgeom_2f1_direct(u, b, bc1, ba1, w, 1, prec);

        acb_rgamma(v, a, prec);
        acb_mul(u, u, v, prec);
        acb_rgamma(v, cb, prec);
        acb_mul(u, u, v, prec);
        acb_neg(v, z);
        acb_neg(cb, b);
        acb_pow(v, v, cb, prec);
        acb_mul(u, u, v, prec);

        acb_sub(t, t, u, prec);

        acb_sin_pi(v, ba, prec);
        acb_div(t, t, v, prec);
        acb_const_pi(v, prec);
        acb_mul(t, t, v, prec);
        acb_set(res, t);

        acb_clear(ca); acb_clear(cb);
        acb_clear(ac1); acb_clear(bc1);
        acb_clear(ab1); acb_clear(ba1);
        acb_clear(w); acb_clear(t);
        acb_clear(u); acb_clear(v);
    }

    acb_clear(ba);
}

