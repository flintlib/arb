/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_elliptic.h"

/* todo: fast handling of y == z (RC function) and other special cases */
/* todo: allow computing cauchy PV */
void
acb_elliptic_rf(acb_t res, const acb_t x, const acb_t y, const acb_t z,
                    int flags, slong prec)
{
    acb_t xx, yy, zz, sx, sy, sz, t;
    acb_t X, Y, Z, E2, E3;
    mag_t err, err2, prev_err;
    slong k, wp, accx, accy, accz;

    if (!acb_is_finite(x) || !acb_is_finite(y) || !acb_is_finite(z))
    {
        acb_indeterminate(res);
        return;
    }

    if (acb_contains_zero(x) + acb_contains_zero(y) + acb_contains_zero(z) > 1)
    {
        acb_indeterminate(res);
        return;
    }

    acb_init(xx); acb_init(yy); acb_init(zz);
    acb_init(sx); acb_init(sy); acb_init(sz);
    acb_init(X); acb_init(Y); acb_init(Z); acb_init(E2); acb_init(E3);
    acb_init(t);
    mag_init(err);
    mag_init(err2);
    mag_init(prev_err);

    acb_set(xx, x);
    acb_set(yy, y);
    acb_set(zz, z);

    wp = prec + 20;

    /* must do at least one iteration */
    for (k = 0; k < prec; k++)
    {
        accx = acb_rel_accuracy_bits(xx);
        accy = acb_rel_accuracy_bits(yy);
        accz = acb_rel_accuracy_bits(zz);

        wp = FLINT_MAX(accx, accy);
        wp = FLINT_MAX(wp, accz);
        wp = FLINT_MAX(wp, 0);
        wp = FLINT_MIN(wp, prec);
        wp += 20;

        acb_sqrt(sx, xx, wp);
        acb_sqrt(sy, yy, wp);
        acb_sqrt(sz, zz, wp);

        acb_add(t, sy, sz, wp);
        acb_mul(t, t, sx, wp);
        acb_addmul(t, sy, sz, wp);

        acb_add(xx, xx, t, wp);
        acb_add(yy, yy, t, wp);
        acb_add(zz, zz, t, wp);

        acb_mul_2exp_si(xx, xx, -2);
        acb_mul_2exp_si(yy, yy, -2);
        acb_mul_2exp_si(zz, zz, -2);

        /* Close enough? Quick estimate based on |x-y|/|x| and |x-z|/|x| */
        /* We also terminate if there is no improvement. */

        acb_sub(t, xx, yy, wp);
        acb_get_mag(err, t);
        acb_sub(t, xx, zz, wp);
        acb_get_mag(err2, t);
        mag_max(err, err, err2);
        acb_get_mag_lower(err2, xx);
        mag_div(err, err, err2);
        mag_pow_ui(err, err, 8);

        if (mag_cmp_2exp_si(err, -prec) < 0 ||
                (k > 2 && mag_cmp(err, prev_err) > 0))
            break;

        mag_set(prev_err, err);
    }

    /* X = 1-x/t, Y = 1-y/t, Z = -X-Y, t = (x+y+z)/3 */
    acb_add(t, xx, yy, prec);
    acb_add(t, t, zz, prec);
    acb_div_ui(t, t, 3, prec);

    acb_div(X, xx, t, prec);
    acb_sub_ui(X, X, 1, prec);
    acb_neg(X, X);

    acb_div(Y, yy, t, prec);
    acb_sub_ui(Y, Y, 1, prec);
    acb_neg(Y, Y);

    acb_add(Z, X, Y, prec);
    acb_neg(Z, Z);

    /* E2 = XY-Z^2, E3 = XYZ */
    acb_mul(E2, X, Y, prec);
    acb_mul(E3, E2, Z, prec);
    acb_submul(E2, Z, Z, prec);

    /*
    Crude bound for the coefficient of
    X^n1 Y^n2 Z^n3 with n1+n2+n3=n: 2*(9/8)^n.
    We truncate before n = 8 (note: for higher precision, we could gain by
    computing more terms). The actual evaluation is done using
    elementary symmetric polynomials, giving

    1 + (-5775*E2**3 + 15015*E2**2*E3 + 10010*E2**2 - 16380*E2*E3 - 24024*E2
         + 6930*E3**2 + 17160*E3)/240240.
    */

    /* Error bound. */
    acb_get_mag(err, X);
    acb_get_mag(err2, Y);
    mag_max(err, err, err2);
    acb_get_mag(err2, Z);
    mag_max(err, err, err2);
    mag_mul_ui(err, err, 9);
    mag_mul_2exp_si(err, err, -3);
    mag_geom_series(err, err, 8);
    mag_mul_2exp_si(err, err, 1);

    /* 1 + (E2*(E2*(-5775*E2 + 15015*E3 + 10010) - 16380*E3 - 24024)
            + E3*(6930*E3 + 17160))/240240 */
    acb_set_ui(sx, 10010);
    acb_addmul_ui(sx, E3, 15015, prec);
    acb_submul_ui(sx, E2, 5775, prec);
    acb_mul(sx, sx, E2, prec);
    acb_submul_ui(sx, E3, 16380, prec);
    acb_sub_ui(sx, sx, 24024, prec);
    acb_mul(sx, sx, E2, prec);
    acb_mul_ui(sy, E3, 6930, prec);
    acb_add_ui(sy, sy, 17160, prec);
    acb_addmul(sx, sy, E3, prec);
    acb_div_ui(sx, sx, 240240, prec);
    acb_add_ui(sx, sx, 1, prec);

    if (acb_is_real(X) && acb_is_real(Y))
        arb_add_error_mag(acb_realref(sx), err);
    else
        acb_add_error_mag(sx, err);

    acb_rsqrt(t, t, prec);
    acb_mul(res, sx, t, prec);

    acb_clear(xx); acb_clear(yy); acb_clear(zz);
    acb_clear(sx); acb_clear(sy); acb_clear(sz);
    acb_clear(X); acb_clear(Y); acb_clear(Z); acb_clear(E2); acb_clear(E3);
    acb_clear(t);
    mag_clear(err);
    mag_clear(err2);
    mag_clear(prev_err);
}

