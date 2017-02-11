/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_elliptic.h"

void
acb_elliptic_rj(acb_t res, const acb_t x, const acb_t y,
            const acb_t z, const acb_t p, int flags, slong prec)
{
    acb_t xx, yy, zz, pp, sx, sy, sz, sp, t, d, delta, S;
    acb_t A, AA, X, Y, Z, P, E2, E3, E4, E5;
    mag_t err, err2;
    slong k;
    int rd;

    if (!acb_is_finite(x) || !acb_is_finite(y) || !acb_is_finite(z) ||
        !acb_is_finite(p))
    {
        acb_indeterminate(res);
        return;
    }

    if ((acb_contains_zero(x) + acb_contains_zero(y) + acb_contains_zero(z) > 1)
        || acb_contains_zero(p))
    {
        acb_indeterminate(res);
        return;
    }

    /* Special case computing R_D(x,y,z) */
    rd = (z == p);

    acb_init(xx); acb_init(yy); acb_init(zz); acb_init(pp);
    acb_init(sx); acb_init(sy); acb_init(sz); acb_init(sp);
    acb_init(S); acb_init(A); acb_init(AA);
    acb_init(X); acb_init(Y); acb_init(Z); acb_init(P);
    acb_init(E2); acb_init(E3); acb_init(E4); acb_init(E5);
    acb_init(t); acb_init(d); acb_init(delta);
    mag_init(err);
    mag_init(err2);

    acb_set(xx, x);
    acb_set(yy, y);
    acb_set(zz, z);
    acb_set(pp, p);
    acb_zero(S);

    if (!rd)
    {
        acb_mul_2exp_si(A, p, 1);
        acb_add(A, A, z, prec);
    }
    else
    {
        acb_mul_ui(A, z, 3, prec);
    }
    acb_add(A, A, x, prec);
    acb_add(A, A, y, prec);
    acb_div_ui(A, A, 5, prec);
    acb_set(AA, A);

    if (!rd)
    {
        acb_sub(delta, p, x, prec);
        acb_sub(t, p, y, prec);
        acb_mul(delta, delta, t, prec);
        acb_sub(t, p, z, prec);
        acb_mul(delta, delta, t, prec);
    }

    /* must do at least one iteration */
    for (k = 0; k < prec; k++)
    {
        acb_sqrt(sx, xx, prec);
        acb_sqrt(sy, yy, prec);
        acb_sqrt(sz, zz, prec);
        if (!rd) acb_sqrt(sp, pp, prec);

        acb_add(t, sy, sz, prec);
        acb_mul(t, t, sx, prec);
        acb_addmul(t, sy, sz, prec);

        acb_add(xx, xx, t, prec);
        acb_add(yy, yy, t, prec);
        acb_add(zz, zz, t, prec);
        if (!rd) acb_add(pp, pp, t, prec);
        acb_add(AA, AA, t, prec);

        acb_mul_2exp_si(xx, xx, -2);
        acb_mul_2exp_si(yy, yy, -2);
        acb_mul_2exp_si(zz, zz, -2);
        if (!rd) acb_mul_2exp_si(pp, pp, -2);
        acb_mul_2exp_si(AA, AA, -2);

        if (!rd)
        {
            /* d = (sp+sx)(sp+sy)(sp+sz) */
            /* e = 4^(-3k) delta / d^2 */
            /* S += 4^(-k) RC(1, 1+e) / d */
            acb_add(d, sp, sx, prec);
            acb_add(t, sp, sy, prec);
            acb_mul(d, d, t, prec);
            acb_add(t, sp, sz, prec);
            acb_mul(d, d, t, prec);

            /* E2 = e */
            acb_mul(E2, d, d, prec);
            acb_div(E2, delta, E2, prec);
            acb_mul_2exp_si(E2, E2, -6 * k);

            acb_elliptic_rc1(E4, E2, prec);
            acb_div(E4, E4, d, prec);
            acb_mul_2exp_si(E4, E4, -2 * k);

            acb_add(S, S, E4, prec);
        }
        else
        {
            acb_mul(t, sz, zz, prec);
            acb_mul_2exp_si(t, t, 2);
            acb_inv(t, t, prec);
            acb_mul_2exp_si(t, t, -2 * k);
            acb_mul_2exp_si(t, t, -1);

            acb_add(S, S, t, prec);
        }

        /* Close enough? */
        acb_sub(t, xx, yy, prec);
        acb_get_mag(err, t);
        acb_sub(t, xx, zz, prec);
        acb_get_mag(err2, t);
        mag_max(err, err, err2);
        if (!rd)
        {
            acb_sub(t, xx, pp, prec);
            acb_get_mag(err2, t);
            mag_max(err, err, err2);
        }
        acb_get_mag_lower(err2, xx);
        mag_div(err, err, err2);

        mag_pow_ui(err, err, 8);

        if (mag_cmp_2exp_si(err, -prec) < 0)
        {
            k++;
            break;
        }
    }

    /* X = (A-x)/(4^k AA) */
    /* Y = (A-y)/(4^k AA) */
    /* Z = (A-z)/(4^k AA) */
    /* P = (-X-Y-Z)/2 */
    acb_mul_2exp_si(t, AA, 2 * k);
    acb_inv(t, t, prec);
    acb_sub(X, A, x, prec);
    acb_mul(X, X, t, prec);
    acb_sub(Y, A, y, prec);
    acb_mul(Y, Y, t, prec);
    acb_sub(Z, A, z, prec);
    acb_mul(Z, Z, t, prec);
    acb_add(P, X, Y, prec);
    acb_add(P, P, Z, prec);
    acb_neg(P, P);
    acb_mul_2exp_si(P, P, -1);

    /* todo: improve for R_D */
    /* E2 = XY + XZ + YZ - 3 P^2 */
    /* E3 = XYZ + 2 E2 P + 4 P^3 */
    /* E4 = (2 XYZ + E2 P + 3 P^3) P */
    /* E5 = XYZP^2 */

    acb_mul(t, P, P, prec); /* t = P^2 */

    acb_mul(E2, X, Y, prec);
    acb_mul(E3, E2, Z, prec);
    acb_mul_2exp_si(E4, E3, 1);
    acb_mul(E5, E3, t, prec);

    acb_add(sx, X, Y, prec);
    acb_addmul(E2, sx, Z, prec);
    acb_submul_ui(E2, t, 3, prec);

    acb_mul(sx, E2, P, prec);
    acb_add(E4, E4, sx, prec);
    acb_mul_2exp_si(sx, sx, 1);
    acb_add(E3, E3, sx, prec);

    acb_mul(t, t, P, prec); /* t = P^3 */
    acb_addmul_ui(E3, t, 4, prec);
    acb_addmul_ui(E4, t, 3, prec);
    acb_mul(E4, E4, P, prec);

    /* Error bound. */
    acb_get_mag(err, X);
    acb_get_mag(err2, Y);
    mag_max(err, err, err2);
    acb_get_mag(err2, Z);
    mag_max(err, err, err2);
    acb_get_mag(err2, P);
    mag_max(err, err, err2);
    mag_mul_ui(err, err, 9);
    mag_mul_2exp_si(err, err, -3);
    mag_geom_series(err, err, 8);
    mag_mul_2exp_si(err, err, 1);

    /*
    1 + (-255255*E2**3 + 675675*E2**2*E3 + 417690*E2**2 - 706860*E2*E3
        + 612612*E2*E4 - 540540*E2*E5 - 875160*E2 + 306306*E3**2 - 540540*E3*E4
        + 680680*E3 - 556920*E4 + 471240*E5)/4084080
    =
    1 + (E2*(E2*(-255255*E2 + 675675*E3 + 417690) - 706860*E3
        + 612612*E4 - 540540*E5 - 875160) + E3*(306306*E3 - 540540*E4
        + 680680) - 556920*E4 + 471240*E5)/4084080
    */
    acb_set_ui(sx, 417690);
    acb_addmul_ui(sx, E3, 675675, prec);
    acb_submul_ui(sx, E2, 255255, prec);
    acb_mul(sx, sx, E2, prec);
    acb_submul_ui(sx, E3, 706860, prec);
    acb_addmul_ui(sx, E4, 612612, prec);
    acb_submul_ui(sx, E5, 540540, prec);
    acb_sub_ui(sx, sx, 875160, prec);
    acb_mul(sx, sx, E2, prec);

    acb_set_ui(sy, 680680);
    acb_submul_ui(sy, E4, 540540, prec);
    acb_addmul_ui(sy, E3, 306306, prec);
    acb_addmul(sx, sy, E3, prec);

    acb_addmul_ui(sx, E5, 471240, prec);
    acb_submul_ui(sx, E4, 556920, prec);
    acb_div_ui(sx, sx, 4084080, prec);

    acb_add_ui(sx, sx, 1, prec);

    if (acb_is_real(X) && acb_is_real(Y) && acb_is_real(Z))
        arb_add_error_mag(acb_realref(sx), err);
    else
        acb_add_error_mag(sx, err);

    acb_rsqrt(t, AA, prec);
    acb_div(t, t, AA, prec);
    acb_mul_2exp_si(t, t, -2 * k);
    acb_mul(t, t, sx, prec);

    acb_addmul_ui(t, S, 6, prec);

    acb_set(res, t);

    acb_clear(xx); acb_clear(yy); acb_clear(zz); acb_clear(pp);
    acb_clear(sx); acb_clear(sy); acb_clear(sz); acb_clear(sp);
    acb_clear(S); acb_clear(A); acb_clear(AA);
    acb_clear(X); acb_clear(Y); acb_clear(Z); acb_clear(P);
    acb_clear(E2); acb_clear(E3); acb_clear(E4); acb_clear(E5);
    acb_clear(t); acb_clear(d); acb_clear(delta);
    mag_clear(err);
    mag_clear(err2);
}

