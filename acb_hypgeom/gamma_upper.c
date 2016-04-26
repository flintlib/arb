/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_hypgeom.h"

void
acb_hypgeom_gamma_upper_asymp(acb_t res, const acb_t s,
        const acb_t z, int regularized, slong prec)
{
    acb_t t, u;

    acb_init(t);
    acb_init(u);

    acb_sub_ui(t, s, 1, prec);
    acb_neg(t, t);

    acb_hypgeom_u_asymp(u, t, t, z, -1, prec);

    if (regularized == 2)
    {
        acb_div(u, u, z, prec);
    }
    else
    {
        acb_neg(t, t);
        acb_pow(t, z, t, prec);
        acb_mul(u, u, t, prec);

        if (regularized == 1)
        {
            acb_rgamma(t, s, prec);
            acb_mul(u, u, t, prec);
        }
    }

    acb_neg(t, z);
    acb_exp(t, t, prec);
    acb_mul(res, t, u, prec);

    acb_clear(t);
    acb_clear(u);
}

void
acb_hypgeom_gamma_upper_1f1a(acb_t res, const acb_t s,
        const acb_t z, int regularized, slong prec)
{
    acb_t a, t, w;
    acb_struct b[2];

    acb_init(a);
    acb_init(b);
    acb_init(b + 1);
    acb_init(t);
    acb_init(w);

    acb_set(a, s);
    acb_add_ui(b, s, 1, prec);
    acb_one(b + 1);
    acb_neg(w, z);

    /* t = 1F1(s, s+1, -z) / s */
    acb_hypgeom_pfq_direct(t, a, 1, b, 2, w, -1, prec);
    acb_div(t, t, s, prec);

    if (regularized == 2)
    {
        acb_neg(a, s);
        acb_pow(a, z, a, prec);
        acb_gamma(b, s, prec);
        acb_mul(b, b, a, prec);
        acb_sub(res, b, t, prec);
    }
    else
    {
        acb_pow(a, z, s, prec);
        acb_mul(t, t, a, prec);

        if (regularized == 1)
        {
            acb_rgamma(a, s, prec);
            acb_mul(t, t, a, prec);
            acb_sub_ui(res, t, 1, prec);
            acb_neg(res, res);
        }
        else
        {
            acb_gamma(a, s, prec);
            acb_sub(res, a, t, prec);
        }
    }

    acb_clear(a);
    acb_clear(b);
    acb_clear(b + 1);
    acb_clear(t);
    acb_clear(w);
}

void
acb_hypgeom_gamma_upper_1f1b(acb_t res, const acb_t s,
        const acb_t z, int regularized, slong prec)
{
    acb_t a, b, t;

    acb_init(a);
    acb_init(b);
    acb_init(t);

    acb_add_ui(b, s, 1, prec);
    acb_hypgeom_pfq_direct(t, NULL, 0, b, 1, z, -1, prec);
    acb_div(t, t, s, prec);

    acb_neg(a, z);
    acb_exp(a, a, prec);
    acb_mul(t, t, a, prec);

    if (regularized == 2)
    {
        acb_neg(a, s);
        acb_pow(a, z, a, prec);
        acb_gamma(b, s, prec);
        acb_mul(b, b, a, prec);
        acb_sub(res, b, t, prec);
    }
    else
    {
        acb_pow(a, z, s, prec);
        acb_mul(t, t, a, prec);

        if (regularized == 1)
        {
            acb_rgamma(a, s, prec);
            acb_mul(t, t, a, prec);
            acb_sub_ui(res, t, 1, prec);
            acb_neg(res, res);
        }
        else
        {
            acb_gamma(a, s, prec);
            acb_sub(res, a, t, prec);
        }
    }

    acb_clear(a);
    acb_clear(b);
    acb_clear(t);
}

/* requires n <= 0 */
void
acb_hypgeom_gamma_upper_singular(acb_t res, slong s, const acb_t z, int regularized, slong prec)
{
    acb_t A, B, C, t, u;
    acb_struct a[2];
    acb_struct b[2];
    arb_t f;
    slong n;

    if (regularized == 1)
    {
        acb_zero(res);
        return;
    }

    n = -s;

    acb_init(A);
    acb_init(B);
    acb_init(C);
    acb_init(t);
    acb_init(u);
    arb_init(f);
    acb_init(a + 0);
    acb_init(a + 1);
    acb_init(b + 0);
    acb_init(b + 1);

    arb_fac_ui(f, n, prec);

    /* (-1)^n (psi(n+1) - log(z)) / n! */
    acb_set_ui(A, n + 1);
    acb_digamma(A, A, prec);

    acb_log(t, z, prec);
    acb_sub(A, A, t, prec);
    acb_div_arb(A, A, f, prec);
    if (n % 2)
        acb_neg(A, A);

    /* (-1)^n z 2F2(1,1;2,2+n;-z) / (n+1)! */
    acb_set_si(a, 1);
    acb_set_si(b, 2);
    acb_set_si(b + 1, 2 + n);
    acb_neg(t, z);
    acb_hypgeom_pfq_direct(B, a, 1, b, 2, t, -1, prec);
    acb_mul(B, B, z, prec);
    arb_mul_ui(f, f, n + 1, prec);
    acb_div_arb(B, B, f, prec);
    if (n % 2)
        acb_neg(B, B);

    /* -sum((-z)^k / ((k-n) k!), k=0...n-1) */
    acb_set_si(a, -n);
    acb_set_si(b, 1 - n);
    acb_set_si(b + 1, 1);
    acb_neg(t, z);

    if (n == 0)
    {
        acb_zero(C);
    }
    else
    {
        acb_hypgeom_pfq_sum(C, u, a, 1, b, 2, t, n, prec);
        acb_div_ui(C, C, n, prec);
    }

    if (regularized == 2)
    {
        acb_pow_si(t, z, n, prec);
        acb_mul(A, A, t, prec);
        acb_mul(B, B, t, prec);
    }
    else
    {
        acb_pow_si(t, z, -n, prec);
        acb_mul(C, C, t, prec);
    }

    acb_add(res, A, B, prec);
    acb_add(res, res, C, prec);

    acb_clear(A);
    acb_clear(B);
    acb_clear(C);
    acb_clear(t);
    acb_clear(u);
    arb_clear(f);
    acb_clear(a + 0);
    acb_clear(a + 1);
    acb_clear(b + 0);
    acb_clear(b + 1);
}

void
acb_hypgeom_gamma_upper(acb_t res, const acb_t s, const acb_t z, int regularized, slong prec)
{
    if (!acb_is_finite(s) || !acb_is_finite(z))
    {
        acb_indeterminate(res);
        return;
    }

    if (acb_is_zero(z))
    {
        if (regularized == 2)
        {
            if (arb_is_negative(acb_realref(s)))
            {
                acb_inv(res, s, prec);
                acb_neg(res, res);
            }
            else
            {
                acb_indeterminate(res);
            }
        }
        else
        {
            if (arb_is_positive(acb_realref(s)))
            {
                if (regularized == 1)
                    acb_one(res);
                else
                    acb_gamma(res, s, prec);
            }
            else
            {
                acb_indeterminate(res);
            }
        }
    }
    else
    {
        slong n = WORD_MAX;

        if (acb_is_int(s))
        {
            if (regularized == 1 && arf_sgn(arb_midref(acb_realref(s))) <= 0)
            {
                acb_zero(res);
                return;
            }

            if (arf_cmpabs_2exp_si(arb_midref(acb_realref(s)), 30) < 0)
            {
                n = arf_get_si(arb_midref(acb_realref(s)), ARF_RND_DOWN);
            }
        }

        if (n >= 1 && n <= 3)
        {
            acb_t t, u;

            acb_init(t);
            acb_init(u);

            if (regularized == 2)
                acb_pow_si(u, z, n, prec);

            if (n == 1)
            {
                acb_neg(res, z);
                acb_exp(res, res, prec);
            }
            else if (n == 2)
            {
                acb_add_ui(t, z, 1, prec);
                acb_neg(res, z);
                acb_exp(res, res, prec);
                acb_mul(res, res, t, prec);
            }
            else if (n == 3)
            {
                acb_add_ui(t, z, 2, prec);
                acb_mul(t, t, z, prec);
                acb_add_ui(t, t, 2, prec);
                acb_neg(res, z);
                acb_exp(res, res, prec);
                acb_mul(res, res, t, prec);
            }

            if (regularized == 2)
                acb_div(res, res, u, prec);
            else if (regularized == 1 && n == 3)
                acb_mul_2exp_si(res, res, -1);

            acb_clear(t);
            acb_clear(u);

            return;
        }

        if (acb_hypgeom_u_use_asymp(z, prec))
        {
            acb_hypgeom_gamma_upper_asymp(res, s, z, regularized, prec);
            return;
        }

        if (n <= 0 && n > -10 * prec)
        {
            acb_hypgeom_gamma_upper_singular(res, n, z, regularized, prec);
            return;
        }

        if (arf_sgn(arb_midref(acb_realref(z))) > 0)
            acb_hypgeom_gamma_upper_1f1b(res, s, z, regularized, prec);
        else
            acb_hypgeom_gamma_upper_1f1a(res, s, z, regularized, prec);
    }
}

