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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "double_extras.h"
#include "fmpcb_mat.h"

long _fmprb_mat_exp_choose_N(const fmpr_t norm, long prec);

void _fmprb_mat_exp_bound(fmpr_t err, const fmpr_t norm, long N);

/* evaluates the truncated Taylor series (assumes no aliasing) */
void
_fmpcb_mat_exp_taylor(fmpcb_mat_t S, const fmpcb_mat_t A, long N, long prec)
{
    if (N == 1)
    {
        fmpcb_mat_one(S);
    }
    else if (N == 2)
    {
        fmpcb_mat_one(S);
        fmpcb_mat_add(S, S, A, prec);
    }
    else if (N == 3)
    {
        fmpcb_mat_t T;
        fmpcb_mat_init(T, fmpcb_mat_nrows(A), fmpcb_mat_nrows(A));
        fmpcb_mat_mul(T, A, A, prec);
        fmpcb_mat_scalar_mul_2exp_si(T, T, -1);
        fmpcb_mat_add(S, A, T, prec);
        fmpcb_mat_one(T);
        fmpcb_mat_add(S, S, T, prec);
        fmpcb_mat_clear(T);
    }
    else
    {
        long i, lo, hi, m, w, dim;
        fmpcb_mat_struct * pows;
        fmpcb_mat_t T, U;
        fmpz_t c, f;

        dim = fmpcb_mat_nrows(A);
        m = n_sqrt(N);
        w = (N + m - 1) / m;

        fmpz_init(c);
        fmpz_init(f);
        pows = flint_malloc(sizeof(fmpcb_mat_t) * (m + 1));
        fmpcb_mat_init(T, dim, dim);
        fmpcb_mat_init(U, dim, dim);

        for (i = 0; i <= m; i++)
        {
            fmpcb_mat_init(pows + i, dim, dim);
            if (i == 0)
                fmpcb_mat_one(pows + i);
            else if (i == 1)
                fmpcb_mat_set(pows + i, A);
            else
                fmpcb_mat_mul(pows + i, pows + i - 1, A, prec);
        }

        fmpcb_mat_zero(S);
        fmpz_one(f);

        for (i = w - 1; i >= 0; i--)
        {
            lo = i * m;
            hi = FLINT_MIN(N - 1, lo + m - 1);

            fmpcb_mat_zero(T);
            fmpz_one(c);

            while (hi >= lo)
            {
                fmpcb_mat_scalar_addmul_fmpz(T, pows + hi - lo, c, prec);
                if (hi != 0)
                    fmpz_mul_ui(c, c, hi);
                hi--;
            }

            fmpcb_mat_mul(U, pows + m, S, prec);
            fmpcb_mat_scalar_mul_fmpz(S, T, f, prec);
            fmpcb_mat_add(S, S, U, prec);
            fmpz_mul(f, f, c);
        }

        fmpcb_mat_scalar_div_fmpz(S, S, f, prec);

        fmpz_clear(c);
        fmpz_clear(f);
        for (i = 0; i <= m; i++)
            fmpcb_mat_clear(pows + i);
        flint_free(pows);
        fmpcb_mat_clear(T);
        fmpcb_mat_clear(U);
    }
}

void
fmpcb_mat_exp(fmpcb_mat_t B, const fmpcb_mat_t A, long prec)
{
    long i, j, dim, wp, N, q, r;
    fmpr_t norm, err;
    fmpcb_mat_t T;

    dim = fmpcb_mat_nrows(A);

    if (dim != fmpcb_mat_ncols(A))
    {
        printf("fmpcb_mat_exp: a square matrix is required!\n");
        abort();
    }

    if (dim == 0)
    {
        return;
    }
    else if (dim == 1)
    {
        fmpcb_exp(fmpcb_mat_entry(B, 0, 0), fmpcb_mat_entry(A, 0, 0), prec);
        return;
    }

    wp = prec + 3 * FLINT_BIT_COUNT(prec);

    fmpr_init(norm);
    fmpr_init(err);
    fmpcb_mat_init(T, dim, dim);

    fmpcb_mat_bound_inf_norm(norm, A, FMPRB_RAD_PREC);

    if (fmpr_is_zero(norm))
    {
        fmpcb_mat_one(B);
    }
    else
    {
        r = fmpr_abs_bound_lt_2exp_si(norm);
        q = pow(wp, 0.25);  /* wanted magnitude */

        if (r > 2 * wp)  /* too big */
            r = 2 * wp;
        else if (r < -q) /* tiny, no need to reduce */
            r = 0;
        else
            r += q;  /* reduce to magnitude 2^(-r) */

        fmpcb_mat_scalar_mul_2exp_si(T, A, -r);
        fmpr_mul_2exp_si(norm, norm, -r);

        N = _fmprb_mat_exp_choose_N(norm, wp);
        _fmprb_mat_exp_bound(err, norm, N);

        _fmpcb_mat_exp_taylor(B, T, N, wp);

        for (i = 0; i < dim; i++)
            for (j = 0; j < dim; j++)
                fmpcb_add_error_fmpr(fmpcb_mat_entry(B, i, j), err);

        for (i = 0; i < r; i++)
        {
            fmpcb_mat_mul(T, B, B, wp);
            fmpcb_mat_swap(T, B);
        }

        for (i = 0; i < dim; i++)
            for (j = 0; j < dim; j++)
                fmpcb_set_round(fmpcb_mat_entry(B, i, j),
                    fmpcb_mat_entry(B, i, j), prec);
    }

    fmpr_clear(norm);
    fmpr_clear(err);
    fmpcb_mat_clear(T);
}

