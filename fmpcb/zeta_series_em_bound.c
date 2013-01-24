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

#include "fmpcb.h"
#include "fmprb_poly.h"
#include "fmpcb_poly.h"

void
fmprb_pow(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_log(t, x, prec);
    fmprb_mul(t, t, y, prec);
    fmprb_exp(z, t, prec);
    fmprb_clear(t);
}

void
bound_I(fmprb_struct * I, const fmprb_t A, const fmprb_t B, const fmprb_t C, long len, long wp)
{
    long k;

    fmprb_t D, Dk, L, T, Bm1;

    fmprb_init(D);
    fmprb_init(Dk);
    fmprb_init(Bm1);
    fmprb_init(T);
    fmprb_init(L);

    fmprb_sub_ui(Bm1, B, 1, wp);
    fmprb_one(L);

    /* T = 1 / (A^Bm1 * Bm1) */
    fmprb_ui_div(T, 1, A, wp);
    fmprb_pow(T, T, Bm1, wp);
    fmprb_div(T, T, Bm1, wp);

    if (len > 1)
    {
        fmprb_log(D, A, wp);
        fmprb_add(D, D, C, wp);
        fmprb_mul(D, D, Bm1, wp);
        fmprb_set(Dk, D);
    }

    for (k = 0; k < len; k++)
    {
        if (k > 0)
        {
            fmprb_mul_ui(L, L, k, wp);
            fmprb_add(L, L, Dk, wp);
            fmprb_mul(Dk, Dk, D, wp);
        }

        fmprb_mul(I + k, L, T, wp);
        fmprb_div(T, T, Bm1, wp);
    }

    fmprb_clear(D);
    fmprb_clear(Dk);
    fmprb_clear(Bm1);
    fmprb_clear(T);
    fmprb_clear(L);
}

/* 0.5*(B/AN)^2 + |B|/AN */
void
bound_C(fmprb_t C, const fmprb_t AN, const fmprb_t B, long wp)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_abs(t, B);
    fmprb_div(t, t, AN, wp);
    fmprb_mul_2exp_si(C, t, -1);
    fmprb_add_ui(C, C, 1, wp);
    fmprb_mul(C, C, t, wp);
    fmprb_clear(t);
}

void
bound_K(fmprb_t C, const fmprb_t AN, const fmprb_t B, const fmprb_t T, long wp)
{
    if (fmprb_is_zero(B) || fmprb_is_zero(T))
    {
        fmprb_one(C);
    }
    else
    {
        fmprb_div(C, B, AN, wp);
        /* TODO: atan is dumb, should also bound by pi/2 */
        fmprb_atan(C, C, wp);
        fmprb_mul(C, C, T, wp);
        if (fmprb_is_nonpositive(C))
            fmprb_one(C);
        else
            fmprb_exp(C, C, wp);
    }
}

/* Absolute value of rising factorial (could speed up once complex gamma is available). */
void
fmpcb_rfac_abs_ubound2(fmpr_t bound, const fmpcb_t s, ulong n, long prec)
{
    fmpr_t term, t;
    ulong k;

    /* M(k) = (a+k)^2 + b^2
       M(0) = a^2 + b^2
       M(k+1) = M(k) + 2*a + (2*k+1)
    */
    fmpr_init(t);
    fmpr_init(term);

    fmpr_one(bound);

    /* M(0) = a^2 + b^2 */
    fmprb_get_abs_ubound_fmpr(t, fmpcb_realref(s), prec);
    fmpr_mul(term, t, t, prec, FMPR_RND_UP);
    fmprb_get_abs_ubound_fmpr(t, fmpcb_imagref(s), prec);
    fmpr_mul(t, t, t, prec, FMPR_RND_UP);
    fmpr_add(term, term, t, prec, FMPR_RND_UP);

    /* we add t = 2*a to each term. note that this can be signed;
       we always want the most positive value */
    fmpr_add(t, fmprb_midref(fmpcb_realref(s)),
        fmprb_radref(fmpcb_realref(s)), prec, FMPR_RND_CEIL);
    fmpr_mul_2exp_si(t, t, 1);

    for (k = 0; k < n; k++)
    {
        fmpr_mul(bound, bound, term, prec, FMPR_RND_UP);
        fmpr_add_ui(term, term, 2 * k + 1, prec, FMPR_RND_UP);
        fmpr_add(term, term, t, prec, FMPR_RND_UP);
    }

    fmpr_sqrt(bound, bound, prec, FMPR_RND_UP);

    fmpr_clear(t);
    fmpr_clear(term);
}


void
bound_rfac(fmprb_struct * F, const fmpcb_t s, ulong n, long len, long wp)
{
    if (len == 1)
    {
        fmpcb_rfac_abs_ubound2(fmprb_midref(F + 0), s, n, wp);
        fmpr_zero(fmprb_radref(F + 0));
    }
    else
    {
        fmprb_struct sx[2];
        fmprb_init(sx + 0);
        fmprb_init(sx + 1);
        fmpcb_abs(sx + 0, s, wp);
        fmprb_one(sx + 1);
        _fmprb_vec_zero(F, len);
        _fmprb_poly_rfac_series_ui(F, sx, 2, n, len, wp);
        fmprb_clear(sx + 0);
        fmprb_clear(sx + 1);
    }
}

void
fmpcb_zeta_series_em_vec_bound(fmprb_struct * bound, const fmpcb_t s, const fmpcb_t a, ulong N, ulong M, long len, long wp)
{
    fmprb_t K, C, AN, S2M;
    fmprb_struct *F, *R;
    long k;

    const fmprb_struct * alpha = fmpcb_realref(a);
    const fmprb_struct * beta  = fmpcb_imagref(a);
    const fmprb_struct * sigma = fmpcb_realref(s);
    const fmprb_struct * tau   = fmpcb_imagref(s);

    fmprb_init(AN);
    fmprb_init(S2M);

    /* require alpha + N > 1, sigma + 2M > 1 */
    fmprb_add_ui(AN, alpha, N - 1, wp);
    fmprb_add_ui(S2M, sigma, 2*M - 1, wp);

    if (!fmprb_is_positive(AN) || !fmprb_is_positive(S2M) || N < 1 || M < 1)
    {
        fmprb_clear(AN);
        fmprb_clear(S2M);

        for (k = 0; k < len; k++)
        {
            fmpr_pos_inf(fmprb_midref(bound + k));
            fmpr_zero(fmprb_radref(bound + k));
        }
        return;
    }

    /* alpha + N, sigma + 2M */
    fmprb_add_ui(AN, AN, 1, wp);
    fmprb_add_ui(S2M, S2M, 1, wp);

    R = _fmprb_vec_init(len);
    F = _fmprb_vec_init(len);

    fmprb_init(K);
    fmprb_init(C);

    /* bound for power integral */
    bound_C(C, AN, beta, wp);
    bound_K(K, AN, beta, tau, wp);
    bound_I(R, AN, S2M, C, len, wp);

    for (k = 0; k < len; k++)
    {
        fmprb_mul(R + k, R + k, K, wp);
        fmprb_div_ui(K, K, k + 1, wp);
    }

    /* bound for rising factorial */
    bound_rfac(F, s, 2*M, len, wp);

    /* product */
    _fmprb_poly_mullow(bound, F, len, R, len, len, wp);

    /* bound for bernoulli polynomials, 4 / (2pi)^(2M) */
    fmprb_const_pi(C, wp);
    fmprb_mul_2exp_si(C, C, 1);
    fmprb_pow_ui(C, C, 2 * M, wp);
    fmprb_ui_div(C, 4, C, wp);
    _fmprb_vec_scalar_mul(bound, bound, len, C, wp);

    fmprb_clear(K);
    fmprb_clear(C);
    fmprb_clear(AN);
    fmprb_clear(S2M);

    _fmprb_vec_clear(R, len);
    _fmprb_vec_clear(F, len);
}

void
fmpcb_zeta_series_em_bound(fmpr_t bound,
        const fmpcb_t s, const fmpcb_t a, long N, long M, long len, long wp)
{
    fmprb_struct * vec = _fmprb_vec_init(len);
    fmpcb_zeta_series_em_vec_bound(vec, s, a, N, M, len, wp);
    _fmprb_vec_get_abs_ubound_fmpr(bound, vec, len, wp);
    _fmprb_vec_clear(vec, len);
}
