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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "fmpz_mat.h"
#include "double_extras.h"
#include "arb_mat.h"
#include "fmpz_mat_extras.h"

#define LOG2_OVER_E 0.25499459743395350926

slong
_fmpz_mat_max_si(const fmpz_mat_t A)
{
    slong i, j, x;

    /* empty max is undefined */
    if (fmpz_mat_is_empty(A)) abort(); /* assert */

    x = fmpz_get_si(fmpz_mat_entry(A, 0, 0));
    for (i = 0; i < fmpz_mat_nrows(A); i++)
    {
        for (j = 0; j < fmpz_mat_ncols(A); j++)
        {
            x = FLINT_MAX(x, fmpz_get_si(fmpz_mat_entry(A, i, j)));
        }
    }

    return x;
}

slong
_arb_mat_exp_choose_N(const mag_t norm, slong prec)
{
    if (mag_is_special(norm) || mag_cmp_2exp_si(norm, 30) > 0 ||
        mag_cmp_2exp_si(norm, -prec) < 0)
    {
        return 1;
    }
    else if (mag_cmp_2exp_si(norm, -300) < 0)
    {
        slong N = -MAG_EXP(norm);
        return (prec + N - 1) / N;
    }
    else
    {
        double c, t;

        c = mag_get_d(norm);
        t = d_lambertw(prec * LOG2_OVER_E / c);
        t = c * exp(t + 1.0);
        return FLINT_MIN((slong) (t + 1.0), 2 * prec);
    }
}

/* evaluates the truncated Taylor series (assumes no aliasing) */
void
_arb_mat_exp_taylor(arb_mat_t S, const arb_mat_t A, slong N, slong prec)
{
    if (N == 1)
    {
        arb_mat_one(S);
    }
    else if (N == 2)
    {
        arb_mat_one(S);
        arb_mat_add(S, S, A, prec);
    }
    else if (N == 3)
    {
        arb_mat_t T;
        arb_mat_init(T, arb_mat_nrows(A), arb_mat_nrows(A));
        arb_mat_sqr(T, A, prec);
        arb_mat_scalar_mul_2exp_si(T, T, -1);
        arb_mat_add(S, A, T, prec);
        arb_mat_one(T);
        arb_mat_add(S, S, T, prec);
        arb_mat_clear(T);
    }
    else
    {
        slong i, lo, hi, m, w, dim;
        arb_mat_struct * pows;
        arb_mat_t T, U;
        fmpz_t c, f;

        dim = arb_mat_nrows(A);
        m = n_sqrt(N);
        w = (N + m - 1) / m;

        fmpz_init(c);
        fmpz_init(f);
        pows = flint_malloc(sizeof(arb_mat_t) * (m + 1));
        arb_mat_init(T, dim, dim);
        arb_mat_init(U, dim, dim);

        for (i = 0; i <= m; i++)
        {
            arb_mat_init(pows + i, dim, dim);
            if (i == 0)
                arb_mat_one(pows + i);
            else if (i == 1)
                arb_mat_set(pows + i, A);
            else
                arb_mat_mul(pows + i, pows + i - 1, A, prec);
        }

        arb_mat_zero(S);
        fmpz_one(f);

        for (i = w - 1; i >= 0; i--)
        {
            lo = i * m;
            hi = FLINT_MIN(N - 1, lo + m - 1);

            arb_mat_zero(T);
            fmpz_one(c);

            while (hi >= lo)
            {
                arb_mat_scalar_addmul_fmpz(T, pows + hi - lo, c, prec);
                if (hi != 0)
                    fmpz_mul_ui(c, c, hi);
                hi--;
            }

            arb_mat_mul(U, pows + m, S, prec);
            arb_mat_scalar_mul_fmpz(S, T, f, prec);
            arb_mat_add(S, S, U, prec);
            fmpz_mul(f, f, c);
        }

        arb_mat_scalar_div_fmpz(S, S, f, prec);

        fmpz_clear(c);
        fmpz_clear(f);
        for (i = 0; i <= m; i++)
            arb_mat_clear(pows + i);
        flint_free(pows);
        arb_mat_clear(T);
        arb_mat_clear(U);
    }
}

void
arb_mat_exp(arb_mat_t B, const arb_mat_t A, slong prec)
{
    slong i, j, dim, wp, N, q, r;
    fmpz_mat_t P;

    if (!arb_mat_is_square(A))
    {
        flint_printf("arb_mat_exp: a square matrix is required!\n");
        abort();
    }

    if (arb_mat_is_empty(A))
        return;

    dim = arb_mat_nrows(A);

    if (dim == 1)
    {
        arb_exp(arb_mat_entry(B, 0, 0), arb_mat_entry(A, 0, 0), prec);
        return;
    }

    fmpz_mat_init(P, dim, dim);
    arb_mat_entrywise_not_is_zero(P, A);

    if (fmpz_mat_is_diagonal(P))
    {
        if (fmpz_mat_is_zero(P))
        {
            arb_mat_one(B);
        }
        else
        {
            if (B != A)
            {
                arb_mat_zero(B);
            }
            for (i = 0; i < dim; i++)
            {
                arb_exp(arb_mat_entry(B, i, i), arb_mat_entry(A, i, i), prec);
            }
        }
    }
    else
    {
        mag_t norm, err;
        arb_mat_t T;
        fmpz_mat_t S;
        int using_structure;

        using_structure = fmpz_mat_count_nonzero(P) < dim * dim;
        if (using_structure)
        {
            fmpz_mat_init(S, dim, dim);
            fmpz_mat_unweighted_all_pairs_longest_walk(S, P);
            fmpz_mat_add_ui_entrywise(S, S, 1);
        }

        mag_init(norm);
        mag_init(err);
        arb_mat_init(T, dim, dim);

        wp = prec + 3 * FLINT_BIT_COUNT(prec);
        arb_mat_bound_inf_norm(norm, A);

        q = pow(wp, 0.25);  /* wanted magnitude */

        if (mag_cmp_2exp_si(norm, 2 * wp) > 0) /* too big */
            r = 2 * wp;
        else if (mag_cmp_2exp_si(norm, -q) < 0) /* tiny, no need to reduce */
            r = 0;
        else
            r = FLINT_MAX(0, q + MAG_EXP(norm)); /* reduce to magnitude 2^(-r) */

        arb_mat_scalar_mul_2exp_si(T, A, -r);
        mag_mul_2exp_si(norm, norm, -r);

        N = _arb_mat_exp_choose_N(norm, wp);

        /* get a second opinion if the matrix is nilpotent */
        if (using_structure && fmpz_mat_is_nonnegative(S))
        {
            N = FLINT_MIN(N, _fmpz_mat_max_si(S));
        }

        mag_exp_tail(err, norm, N);

        _arb_mat_exp_taylor(B, T, N, wp);

        if (using_structure)
        {
            for (i = 0; i < dim; i++)
            {
                for (j = 0; j < dim; j++)
                {
                    if (fmpz_sgn(fmpz_mat_entry(S, i, j)) < 0 ||
                        fmpz_cmp_si(fmpz_mat_entry(S, i, j), N) > 0)
                    {
                        arb_add_error_mag(arb_mat_entry(B, i, j), err);
                    }
                }
            }
            fmpz_mat_clear(S);
        }
        else
        {
            for (i = 0; i < dim; i++)
                for (j = 0; j < dim; j++)
                    arb_add_error_mag(arb_mat_entry(B, i, j), err);
        }

        for (i = 0; i < r; i++)
        {
            arb_mat_sqr(T, B, wp);
            arb_mat_swap(T, B);
        }

        for (i = 0; i < dim; i++)
            for (j = 0; j < dim; j++)
                arb_set_round(arb_mat_entry(B, i, j),
                    arb_mat_entry(B, i, j), prec);

        mag_clear(norm);
        mag_clear(err);
        arb_mat_clear(T);
    }

    fmpz_mat_clear(P);
}
