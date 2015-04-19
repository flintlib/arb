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
_acb_poly_reciprocal_majorant(arb_ptr res, acb_srcptr vec, long len, long prec)
{
    long i;

    for (i = 0; i < len; i++)
    {
        if (i == 0)
        {
            acb_get_abs_lbound_arf(arb_midref(res + i), vec + i, prec);
            mag_zero(arb_radref(res + i));
        }
        else
        {
            acb_get_abs_ubound_arf(arb_midref(res + i), vec + i, prec);
            arf_neg(arb_midref(res + i), arb_midref(res + i));
            mag_zero(arb_radref(res + i));
        }
    }
}

void
acb_poly_reciprocal_majorant(arb_poly_t res, const acb_poly_t poly, long prec)
{
    arb_poly_fit_length(res, poly->length);
    _acb_poly_reciprocal_majorant(res->coeffs, poly->coeffs, poly->length, prec);
    _arb_poly_set_length(res, poly->length);
}


/* F = 1 + U + U^2 + ... = 1/(1-U) assuming that U[0] is positive;
   indeterminate if not convergent */
static void
arb_poly_geometric_sum(arb_poly_t F, const arb_poly_t U, long len, long prec)
{
    if (U->length == 0)
    {
        arb_poly_one(F);
        return;
    }

    arb_poly_add_si(F, U, -1, prec);
    arb_poly_neg(F, F);

    if (F->length > 0 && arb_is_positive(F->coeffs))
    {
        arb_poly_inv_series(F, F, len, prec);
    }
    else
    {
        arb_poly_fit_length(F, len);
        _arb_vec_indeterminate(F->coeffs, len);
        _arb_poly_set_length(F,  len);
    }
}

/* F = 1 + U + U^2 + U^3 + ... = 1/(1-U)

   U = product of (1 + |A-B|/(|B[0] - |B[1:]|)
       product of (1 / (|B[0] - |B[1:]|))
       * |Z|
*/
void
acb_hypgeom_pfq_series_bound_factor(arb_poly_t F,
    const acb_poly_struct * a, long p,
    const acb_poly_struct * b, long q, 
    const acb_poly_t z,
    long n, long len, long prec)
{
    long i;

    arb_poly_t T, U, V;
    acb_poly_t BN, AB;

    /* not convergent */
    if (p > q)
    {
        arb_poly_fit_length(F, len);
        _arb_vec_indeterminate(F->coeffs, len);
        _arb_poly_set_length(F, len);
        return;
    }

    arb_poly_init(T);
    arb_poly_init(U);
    arb_poly_init(V);

    acb_poly_init(BN);
    acb_poly_init(AB);

    acb_poly_majorant(U, z, prec);

    for (i = 0; i < q; i++)
    {
        acb_poly_add_si(BN, b + i, n, prec);

        if (acb_poly_length(BN) != 0 &&
                arb_is_positive(acb_realref(BN->coeffs)))
        {
            if (i < p)
            {
                /* 1 + |a-b|/reciprocal_majorant(b + n) */
                acb_poly_sub(AB, a + i, b + i, prec);
                acb_poly_majorant(T, AB, prec);
                acb_poly_reciprocal_majorant(V, BN, prec);
                arb_poly_div_series(T, T, V, len, prec);
                arb_poly_add_si(T, T, 1, prec);
                arb_poly_mullow(U, U, T, len, prec);
            }
            else
            {
                acb_poly_reciprocal_majorant(T, BN, prec);
                arb_poly_div_series(U, U, T, len, prec);
            }
        }
        else
        {
            arb_poly_fit_length(U, len);
            _arb_vec_indeterminate(U->coeffs, len);
            _arb_poly_set_length(U,  len);
            break;
        }
    }

    /* F = 1/(1-U) */
    arb_poly_geometric_sum(F, U, len, prec);

    arb_poly_clear(T);
    arb_poly_clear(U);
    arb_poly_clear(V);

    acb_poly_clear(BN);
    acb_poly_clear(AB);
}

void
acb_hypgeom_pfq_series_sum_forward(acb_poly_t s, acb_poly_t t,
    const acb_poly_struct * a, long p,
    const acb_poly_struct * b, long q,
    const acb_poly_t z, int regularized,
    long n, long len, long prec)
{
    acb_poly_t u, v;
    acb_poly_t tmp;
    long k, i;

    acb_poly_init(u);
    acb_poly_init(v);
    acb_poly_init(tmp);

    if (!regularized)
    {
        acb_poly_zero(s);
        acb_poly_one(t);

        for (k = 0; k < n && acb_poly_length(t) != 0; k++)
        {
            acb_poly_add(s, s, t, prec);

            if (p > 0)
            {
                acb_poly_add_si(u, a, k, prec);

                for (i = 1; i < p; i++)
                {
                    acb_poly_add_si(v, a + i, k, prec);
                    acb_poly_mullow(u, u, v, len, prec);
                }

                acb_poly_mullow(t, t, u, len, prec);
            }

            if (q > 0)
            {
                acb_poly_add_si(u, b, k, prec);

                for (i = 1; i < q; i++)
                {
                    acb_poly_add_si(v, b + i, k, prec);
                    acb_poly_mullow(u, u, v, len, prec);
                }

                acb_poly_div_series(t, t, u, len, prec);
            }

            acb_poly_mullow(t, t, z, len, prec);
        }
    }
    else
    {
        acb_poly_zero(s);

        for (i = 0; i < q; i++)
        {
            if (i == 0)
            {
                acb_poly_rgamma_series(t, b + i, len, prec);
            }
            else
            {
                acb_poly_rgamma_series(u, b + i, len, prec);
                acb_poly_mullow(tmp, t, u, len, prec);
                acb_poly_swap(tmp, t);
            }
        }

        for (k = 0; k < n; k++)
        {
            acb_poly_add(s, s, t, prec);

            if (p > 0)
            {
                acb_poly_add_si(u, a, k, prec);

                for (i = 1; i < p; i++)
                {
                    acb_poly_add_si(v, a + i, k, prec);
                    acb_poly_mullow(tmp, u, v, len, prec);
                    acb_poly_swap(tmp, u);
                }

                acb_poly_mullow(tmp, t, u, len, prec);
                acb_poly_swap(tmp, t);
            }

            if (q > 0)
            {
                acb_poly_add_si(u, b, k, prec);

                for (i = 1; i < q; i++)
                {
                    acb_poly_add_si(v, b + i, k, prec);
                    acb_poly_mullow(tmp, u, v, len, prec);
                    acb_poly_swap(tmp, u);
                }

                if (acb_poly_length(u) > 0 && !acb_contains_zero(u->coeffs))
                {
                    acb_poly_div_series(tmp, t, u, len, prec);
                    acb_poly_mullow(t, tmp, z, len, prec);
                }
                else
                {
                    /* compute term from scratch */
                    acb_poly_one(t);

                    for (i = 0; i < p; i++)
                    {
                        acb_poly_rising_ui_series(v, a + i, k + 1, len, prec);
                        acb_poly_mullow(t, t, v, len, prec);
                    }

                    for (i = 0; i < q; i++)
                    {
                        acb_poly_add_si(v, b + i, k + 1, prec);
                        acb_poly_rgamma_series(v, v, len, prec);
                        acb_poly_mullow(t, t, v, len, prec);
                    }

                    acb_poly_pow_ui_trunc_binexp(v, z, k + 1, len, prec);
                    acb_poly_mullow(t, t, v, len, prec);
                }
            }
            else
            {
                acb_poly_mullow(tmp, t, z, len, prec);
                acb_poly_swap(tmp, t);
            }
        }
    }

    acb_poly_clear(u);
    acb_poly_clear(v);
    acb_poly_clear(tmp);
}

void
acb_hypgeom_pfq_series_direct(acb_poly_t res,
    const acb_poly_struct * a, long p,
    const acb_poly_struct * b, long q,
    const acb_poly_t z, int regularized,
    long n, long len, long prec)
{
    acb_poly_t s, t, err;
    arb_poly_t C, T;
    long i;
    int is_real;

    /* default algorithm to choose number of terms */
    if (n < 0)
    {
        n = acb_hypgeom_pfq_series_choose_n(a, p, b, q, z, len, prec);
    }

    acb_poly_init(s);
    acb_poly_init(t);
    acb_poly_init(err);
    arb_poly_init(C);
    arb_poly_init(T);

    acb_hypgeom_pfq_series_sum_forward(s, t, a, p, b, q, z, regularized, n, len, prec);

    if (acb_poly_length(t) != 0)
    {
        is_real = acb_poly_is_real(z);
        for (i = 0; i < p; i++)
            is_real = is_real && acb_poly_is_real(a + i);
        for (i = 0; i < q; i++)
            is_real = is_real && acb_poly_is_real(b + i);

        acb_poly_majorant(T, t, MAG_BITS);
        acb_hypgeom_pfq_series_bound_factor(C, a, p, b, q, z, n, len, MAG_BITS);

        arb_poly_mullow(T, T, C, len, MAG_BITS);

        /* create polynomial of errors */
        acb_poly_fit_length(err, len);

        for (i = 0; i < FLINT_MIN(len, T->length); i++)
        {
            arb_add_error(acb_realref(err->coeffs + i), T->coeffs + i);
            if (!is_real)
                arb_add_error(acb_imagref(err->coeffs + i), T->coeffs + i);
        }

        _acb_poly_set_length(err, len);
        _acb_poly_normalise(err);

        acb_poly_add(s, s, err, prec);
    }

    acb_poly_set(res, s);

    acb_poly_clear(s);
    acb_poly_clear(t);
    acb_poly_clear(err);
    arb_poly_clear(C);
    arb_poly_clear(T);
}

