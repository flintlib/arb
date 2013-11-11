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

#include "double_extras.h"
#include "fmprb_mat.h"

#define LOG2_OVER_E 0.25499459743395350926

void fmpr_gamma_ui_lbound(fmpr_t x, ulong n, long prec);

long
_fmprb_mat_exp_choose_N(const fmpr_t norm, long prec)
{
    if (fmpr_is_special(norm) || fmpr_cmp_2exp_si(norm, 30) > 0)
    {
        return 1;
    }
    else if (fmpr_cmp_2exp_si(norm, -300) < 0)
    {
        long N = -fmpr_abs_bound_lt_2exp_si(norm);

        return (prec + N - 1) / N;
    }
    else
    {
        double c, t;

        c = fmpr_get_d(norm, FMPR_RND_UP);

        t = d_lambertw(prec * LOG2_OVER_E / c);
        t = c * exp(t + 1.0);

        return FLINT_MIN((long) (t + 1.0), 2 * prec);
    }
}

void
_fmprb_mat_exp_bound(fmpr_t err, const fmpr_t norm, long N)
{
    fmpr_t t, u;

    fmpr_init(t);
    fmpr_init(u);

    fmpr_set_si_2exp_si(t, N, -1);

    /* bound by geometric series when N >= 2*c  <=> N/2 >= c */
    if (N > 0 && fmpr_cmp(t, norm) >= 0)
    {
        /* 2 c^N / N! */
        fmpr_pow_sloppy_ui(t, norm, N, FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_gamma_ui_lbound(u, N + 1, FMPRB_RAD_PREC);
        fmpr_div(err, t, u, FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_mul_2exp_si(err, err, 1);
    }
    else
    {
        fmpr_pos_inf(err);
    }

    fmpr_clear(t);
    fmpr_clear(u);
}


/* evaluates the truncated Taylor series (assumes no aliasing) */
void
_fmprb_mat_exp_taylor(fmprb_mat_t S, const fmprb_mat_t A, long N, long prec)
{
    if (N == 1)
    {
        fmprb_mat_one(S);
    }
    else if (N == 2)
    {
        fmprb_mat_one(S);
        fmprb_mat_add(S, S, A, prec);
    }
    else if (N == 3)
    {
        fmprb_mat_t T;
        fmprb_mat_init(T, fmprb_mat_nrows(A), fmprb_mat_nrows(A));
        fmprb_mat_mul(T, A, A, prec);
        fmprb_mat_scalar_mul_2exp_si(T, T, -1);
        fmprb_mat_add(S, A, T, prec);
        fmprb_mat_one(T);
        fmprb_mat_add(S, S, T, prec);
        fmprb_mat_clear(T);
    }
    else
    {
        long i, lo, hi, m, w, dim;
        fmprb_mat_struct * pows;
        fmprb_mat_t T, U;
        fmpz_t c, f;

        dim = fmprb_mat_nrows(A);
        m = n_sqrt(N);
        w = (N + m - 1) / m;

        fmpz_init(c);
        fmpz_init(f);
        pows = flint_malloc(sizeof(fmprb_mat_t) * (m + 1));
        fmprb_mat_init(T, dim, dim);
        fmprb_mat_init(U, dim, dim);

        for (i = 0; i <= m; i++)
        {
            fmprb_mat_init(pows + i, dim, dim);
            if (i == 0)
                fmprb_mat_one(pows + i);
            else if (i == 1)
                fmprb_mat_set(pows + i, A);
            else
                fmprb_mat_mul(pows + i, pows + i - 1, A, prec);
        }

        fmprb_mat_zero(S);
        fmpz_one(f);

        for (i = w - 1; i >= 0; i--)
        {
            lo = i * m;
            hi = FLINT_MIN(N - 1, lo + m - 1);

            fmprb_mat_zero(T);
            fmpz_one(c);

            while (hi >= lo)
            {
                fmprb_mat_scalar_addmul_fmpz(T, pows + hi - lo, c, prec);
                if (hi != 0)
                    fmpz_mul_ui(c, c, hi);
                hi--;
            }

            fmprb_mat_mul(U, pows + m, S, prec);
            fmprb_mat_scalar_mul_fmpz(S, T, f, prec);
            fmprb_mat_add(S, S, U, prec);
            fmpz_mul(f, f, c);
        }

        fmprb_mat_scalar_div_fmpz(S, S, f, prec);

        fmpz_clear(c);
        fmpz_clear(f);
        for (i = 0; i <= m; i++)
            fmprb_mat_clear(pows + i);
        flint_free(pows);
        fmprb_mat_clear(T);
        fmprb_mat_clear(U);
    }
}

void
fmprb_mat_exp(fmprb_mat_t B, const fmprb_mat_t A, long prec)
{
    long i, j, dim, wp, N, q, r;
    fmpr_t norm, err;
    fmprb_mat_t T;

    dim = fmprb_mat_nrows(A);

    if (dim != fmprb_mat_ncols(A))
    {
        printf("fmprb_mat_exp: a square matrix is required!\n");
        abort();
    }

    if (dim == 0)
    {
        return;
    }
    else if (dim == 1)
    {
        fmprb_exp(fmprb_mat_entry(B, 0, 0), fmprb_mat_entry(A, 0, 0), prec);
        return;
    }

    wp = prec + 3 * FLINT_BIT_COUNT(prec);

    fmpr_init(norm);
    fmpr_init(err);
    fmprb_mat_init(T, dim, dim);

    fmprb_mat_bound_inf_norm(norm, A, FMPRB_RAD_PREC);

    if (fmpr_is_zero(norm))
    {
        fmprb_mat_one(B);
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

        fmprb_mat_scalar_mul_2exp_si(T, A, -r);
        fmpr_mul_2exp_si(norm, norm, -r);

        N = _fmprb_mat_exp_choose_N(norm, wp);
        _fmprb_mat_exp_bound(err, norm, N);

        _fmprb_mat_exp_taylor(B, T, N, wp);

        for (i = 0; i < dim; i++)
            for (j = 0; j < dim; j++)
                fmprb_add_error_fmpr(fmprb_mat_entry(B, i, j), err);

        for (i = 0; i < r; i++)
        {
            fmprb_mat_mul(T, B, B, wp);
            fmprb_mat_swap(T, B);
        }

        for (i = 0; i < dim; i++)
            for (j = 0; j < dim; j++)
                fmprb_set_round(fmprb_mat_entry(B, i, j),
                    fmprb_mat_entry(B, i, j), prec);
    }

    fmpr_clear(norm);
    fmpr_clear(err);
    fmprb_mat_clear(T);
}

