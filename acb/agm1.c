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

#include "acb.h"
#include "acb_poly.h"

/* Checks that |arg(z)| <= 3 pi / 4 */
static int
acb_check_arg(const acb_t z)
{
    mag_t re, im;
    int res;

    if (!arb_contains_negative(acb_realref(z)))
        return 1;

    mag_init(re);
    mag_init(im);

    arb_get_mag(re, acb_realref(z));
    arb_get_mag_lower(im, acb_imagref(z));

    res = mag_cmp(re, im) < 0;

    mag_clear(re);
    mag_clear(im);

    return res;
}

static void
sqrtmul(acb_t c, const acb_t a, const acb_t b, long prec)
{
    if (arb_is_positive(acb_realref(a)) &&
        arb_is_positive(acb_realref(b)))
    {
        acb_mul(c, a, b, prec);
        acb_sqrt(c, c, prec);
    }
    else if (arb_is_nonnegative(acb_imagref(a)) &&
             arb_is_nonnegative(acb_imagref(b)))
    {
        acb_mul(c, a, b, prec);
        acb_neg(c, c);
        acb_sqrt(c, c, prec);
        acb_mul_onei(c, c);
    }
    else if (arb_is_nonpositive(acb_imagref(a)) &&
             arb_is_nonpositive(acb_imagref(b)))
    {
        acb_mul(c, a, b, prec);
        acb_neg(c, c);
        acb_sqrt(c, c, prec);
        acb_mul_onei(c, c);
        acb_neg(c, c);
    }
    else
    {
        acb_t d;
        acb_init(d);
        acb_sqrt(c, a, prec);
        acb_sqrt(d, b, prec);
        acb_mul(c, c, d, prec);
        acb_clear(d);
    }
}

void
acb_agm1_basecase(acb_t res, const acb_t z, long prec)
{
    acb_t a, b, t, u;
    mag_t err;
    int isreal;

    if (acb_is_zero(z))
    {
        acb_zero(res);
        return;
    }

    if (acb_is_one(z))
    {
        acb_one(res);
        return;
    }

    if (!acb_check_arg(z))
    {
        mag_t one;
        mag_init(one);
        mag_init(err);
        mag_one(one);
        acb_get_mag(err, z);
        mag_max(err, err, one);
        acb_zero(res);
        acb_add_error_mag(res, err);
        mag_clear(err);
        mag_clear(one);
        return;
    }

    isreal = acb_is_real(z) && arb_is_nonnegative(acb_realref(z));

    acb_init(a);
    acb_init(b);
    acb_init(t);
    acb_init(u);
    mag_init(err);

    acb_one(a);
    acb_set_round(b, z, prec);

    while (!acb_overlaps(a, b))
    {
        acb_add(t, a, b, prec);
        acb_mul_2exp_si(t, t, -1);

        sqrtmul(u, a, b, prec);

        acb_swap(t, a);
        acb_swap(u, b);
    }

    /* Dupont's thesis, p. 87:
       |M(z) - a_n| <= |a_n - b_n| */
    acb_sub(t, a, b, prec);
    acb_get_mag(err, t);

    if (isreal)
        arb_add_error_mag(acb_realref(a), err);
    else
        acb_add_error_mag(a, err);

    acb_set(res, a);

    acb_clear(a);
    acb_clear(b);
    acb_clear(t);
    acb_clear(u);
    mag_clear(err);
}

/*
    Computes (M(z), M'(z)) using a finite difference.
    Assumes z exact, |arg(z)| <= 3 pi / 4.
*/
void
acb_agm1_deriv_diff(acb_t Mz, acb_t Mzp, const acb_t z, long prec)
{
    mag_t err, t;
    fmpz_t rexp, hexp;
    long wp;
    int isreal;

    if (!acb_is_exact(z) || !acb_is_finite(z) ||
          acb_is_zero(z) || !acb_check_arg(z))
    {
        acb_indeterminate(Mz);
        acb_indeterminate(Mzp);
        return;
    }

    isreal = acb_is_real(z) && arb_is_nonnegative(acb_realref(z));

    /*
       |M^(k)(z) / k!| <= C * D^k  where
       C = max(1, |z| + r),
       D = 1/r,
       and 0 < r < |z|

        M(z+h) - M(z)
       |------------- - M'(z)|  <=  C D^2 h / (1 - D h)
               h

       h D < 1.
    */

    fmpz_init(hexp);
    fmpz_init(rexp);
    mag_init(err);
    mag_init(t);

    /* choose r = 2^rexp such that r < |z| */
    acb_get_mag_lower(t, z);
    fmpz_sub_ui(rexp, MAG_EXPREF(t), 2);

    /* Choose h = 2^hexp with hexp = rexp - (prec + 5).
       D = 2^-rexp
       C D^2 h / (1 - D h) <= C * 2^(-5-prec-rexp+1)
    */

    /* err = C = max(1, |z| + r) */
    acb_get_mag(err, z);
    mag_one(t);
    mag_mul_2exp_fmpz(t, t, rexp);
    mag_add(err, err, t);
    mag_one(t);
    mag_max(err, err, t);

    /* multiply by 2^(-5-prec-rexp+1) (use hexp as temp) */
    fmpz_set_si(hexp, 1 - 5 - prec);
    fmpz_sub(hexp, hexp, rexp);
    mag_mul_2exp_fmpz(err, err, hexp);

    /* choose h = 2^hexp */
    fmpz_sub_ui(hexp, rexp, prec + 5);

    /* compute finite difference */
    wp = 2 * prec + 10;
    acb_agm1_basecase(Mz, z, wp);
    acb_one(Mzp);
    acb_mul_2exp_fmpz(Mzp, Mzp, hexp);
    acb_add(Mzp, Mzp, z, wp);
    acb_agm1_basecase(Mzp, Mzp, wp);
    acb_sub(Mzp, Mzp, Mz, prec);
    fmpz_neg(hexp, hexp);
    acb_mul_2exp_fmpz(Mzp, Mzp, hexp);

    if (isreal)
        arb_add_error_mag(acb_realref(Mzp), err);
    else
        acb_add_error_mag(Mzp, err);

    acb_set_round(Mz, Mz, prec);

    fmpz_clear(hexp);
    fmpz_clear(rexp);
    mag_clear(err);
    mag_clear(t);
}

/*
For input z + eps

First derivative bound: max(1, |z|+|eps|+r) / r
Second derivative bound: 2 max(1, |z|+|eps|+r) / r^2

This is assuming that the circle at z with radius |eps| + r
does not cross the negative half axis, which we check.
*/

void
acb_agm1_deriv_right(acb_t Mz, acb_t Mzp, const acb_t z, long prec)
{
    if (acb_is_exact(z))
    {
        acb_agm1_deriv_diff(Mz, Mzp, z, prec);
    }
    else
    {
        if (!acb_is_finite(z) || !acb_check_arg(z))
        {
            acb_indeterminate(Mz);
            acb_indeterminate(Mzp);
        }
        else
        {
            acb_t t;
            mag_t r, eps, err, one;
            int isreal;

            acb_init(t);
            mag_init(r);
            mag_init(err);
            mag_init(one);
            mag_init(eps);

            isreal = acb_is_real(z) && arb_is_nonnegative(acb_realref(z));

            mag_hypot(eps, arb_radref(acb_realref(z)), arb_radref(acb_imagref(z)));

            /* choose r avoiding overlap with negative half axis */
            if (arf_sgn(arb_midref(acb_realref(z))) < 0)
                arb_get_mag_lower(r, acb_imagref(z));
            else
                acb_get_mag_lower(r, z);

            mag_mul_2exp_si(r, r, -1);

            if (mag_is_zero(r))
            {
                acb_indeterminate(Mz);
                acb_indeterminate(Mzp);
            }
            else
            {
                acb_set(t, z);
                mag_zero(arb_radref(acb_realref(t)));
                mag_zero(arb_radref(acb_imagref(t)));

                acb_get_mag(err, z);
                mag_add(err, err, r);
                mag_add(err, err, eps);
                mag_one(one);
                mag_max(err, err, one);
                mag_mul(err, err, eps);

                acb_agm1_deriv_diff(Mz, Mzp, t, prec);

                mag_div(err, err, r);

                if (isreal)
                    arb_add_error_mag(acb_realref(Mz), err);
                else
                    acb_add_error_mag(Mz, err);

                mag_div(err, err, r);
                mag_mul_2exp_si(err, err, 1);

                if (isreal)
                    arb_add_error_mag(acb_realref(Mzp), err);
                else
                    acb_add_error_mag(Mzp, err);
            }

            acb_clear(t);
            mag_clear(r);
            mag_clear(err);
            mag_clear(one);
            mag_clear(eps);
        }
    }
}

void
acb_agm1(acb_t m, const acb_t z, long prec)
{
    if (arf_sgn(arb_midref(acb_realref(z))) >= 0)
    {
        acb_agm1_basecase(m, z, prec);
    }
    else
    {
        /* use M(z) = (z+1)/2 * M(2 sqrt(z) / (z+1)) */
        acb_t t;
        acb_init(t);
        acb_add_ui(t, z, 1, prec);
        acb_sqrt(m, z, prec);
        acb_div(m, m, t, prec);
        acb_mul_2exp_si(m, m, 1);
        acb_agm1_basecase(m, m, prec);
        acb_mul(m, m, t, prec);
        acb_mul_2exp_si(m, m, -1);
        acb_clear(t);
    }
}

void
acb_agm1_deriv(acb_t Mz, acb_t Mzp, const acb_t z, long prec)
{
    /*
       u = 2 sqrt(z) / (1+z)

       Mz = (1+z) M(u) / 2
       Mzp = [M(u) - (z-1) M'(u) / ((1+z) sqrt(z))] / 2
    */

    if (arf_sgn(arb_midref(acb_realref(z))) >= 0)
    {
        acb_agm1_deriv_right(Mz, Mzp, z, prec);
    }
    else
    {
        acb_t t, u, zp1, zm1;

        acb_init(t);
        acb_init(u);
        acb_init(zp1);
        acb_init(zm1);

        acb_sqrt(t, z, prec);
        acb_add_ui(zp1, z, 1, prec);
        acb_sub_ui(zm1, z, 1, prec);

        acb_div(u, t, zp1, prec);
        acb_mul_2exp_si(u, u, 1);

        acb_agm1_deriv_right(Mz, Mzp, u, prec);

        acb_mul(Mzp, Mzp, zm1, prec);
        acb_mul(t, t, zp1, prec);
        acb_div(Mzp, Mzp, t, prec);
        acb_sub(Mzp, Mz, Mzp, prec);
        acb_mul_2exp_si(Mzp, Mzp, -1);

        acb_mul(Mz, Mz, zp1, prec);
        acb_mul_2exp_si(Mz, Mz, -1);

        acb_clear(t);
        acb_clear(u);
        acb_clear(zp1);
        acb_clear(zm1);
    }
}

void
acb_agm1_cpx(acb_ptr m, const acb_t z, long len, long prec)
{
    if (len < 1)
        return;

    if (len == 1)
    {
        acb_agm1(m, z, prec);
        return;
    }

    if (len == 2)
    {
        acb_agm1_deriv(m, m + 1, z, prec);
        return;
    }

    if (len >= 3)
    {
        acb_t t, u, v;
        acb_ptr w;
        long k, n;

        acb_init(t);
        acb_init(u);
        acb_init(v);
        w = _acb_vec_init(len);

        acb_agm1_deriv(w, w + 1, z, prec);

        /* invert series */
        acb_inv(w, w, prec);
        acb_mul(t, w, w, prec);
        acb_mul(w + 1, w + 1, t, prec);
        acb_neg(w + 1, w + 1);

        if (acb_is_one(z))
        {
            for (k = 2; k < len; k++)
            {
                n = k - 2;

                acb_mul_ui(w + k, w + n + 0, (n+1)*(n+1), prec);
                acb_addmul_ui(w + k, w + n + 1, 7+3*n*(3+n), prec);
                acb_div_ui(w + k, w + k, 2*(n+2)*(n+2), prec);
                acb_neg(w + k, w + k);
            }
        }
        else
        {
            /* t = 3z^2 - 1 */
            /* u = -1 / (z^3 - z) */
            acb_mul(t, z, z, prec);
            acb_mul(u, t, z, prec);
            acb_mul_ui(t, t, 3, prec);
            acb_sub_ui(t, t, 1, prec);
            acb_sub(u, u, z, prec);
            acb_inv(u, u, prec);
            acb_neg(u, u);

            /* use differential equation for second derivative */
            acb_mul(w + 2, z, w + 0, prec);
            acb_addmul(w + 2, t, w + 1, prec);
            acb_mul(w + 2, w + 2, u, prec);
            acb_mul_2exp_si(w + 2, w + 2, -1);

            /* recurrence */
            for (k = 3; k < len; k++)
            {
                n = k - 3;
                acb_mul_ui(w + k, w + n + 0, (n+1)*(n+1), prec);
                acb_mul(v, w + n + 1, z, prec);
                acb_addmul_ui(w + k, v, 7+3*n*(3+n), prec);
                acb_mul(v, w + n + 2, t, prec);
                acb_addmul_ui(w + k, v, (n+2)*(n+2), prec);
                acb_mul(w + k, w + k, u, prec);
                acb_div_ui(w + k, w + k, (n+2)*(n+3), prec);
            }
        }

        /* invert series */
        _acb_poly_inv_series(m, w, len, len, prec);

        acb_clear(t);
        acb_clear(u);
        acb_clear(v);
        _acb_vec_clear(w, len);
    }
}

