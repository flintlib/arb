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

    Copyright (C) 2012, 2013 Fredrik Johansson

******************************************************************************/

#include "fmprb.h"
#include "fmprb_poly.h"
#include "gamma.h"

static void
gamma_sum_bsplit(fmprb_t P, fmprb_t Q,
            fmprb_struct * B, fmprb_struct * T,
            const fmpz_t p, const fmpz_t q,
            long a, long b, long N, long L, long wp)
{
    if (b - a == 1)
    {
        /* P = T = -N, Q = a, B = (q*a + p) + q*x */
        fmprb_set_si(P, (a == 0) ? 1 : -N);
        fmprb_set_si(Q, (a == 0) ? 1 : a);
        fmprb_set(T, P);
        fmprb_set_fmpz(B, q);
        fmprb_mul_ui(B, B, a, wp);
        fmprb_add_fmpz(B, B, p, wp);
        if (L > 1)
            fmprb_set_fmpz(B + 1, q);
    }
    else
    {
        long m, lenbl, lentl, lenbr, lentr, lent, lenb, alloc;
        fmprb_t PL, QL, PR, QR;
        fmprb_struct *BL, *TL, *BR, *TR, *TT;

        fmprb_init(PL);
        fmprb_init(QL);
        fmprb_init(PR);
        fmprb_init(QR);

        m = a + (b - a) / 2;

        lenbl = FLINT_MIN(m - a + 1, L);
        lenbr = FLINT_MIN(b - m + 1, L);
        lentl = FLINT_MIN(m - a, L);
        lentr = FLINT_MIN(b - m, L);
        lenb = FLINT_MIN(b - a + 1, L);
        lent = FLINT_MIN(b - a, L);
        alloc = lenbl + lenbr + lentl + lentr + lent;

        BL = _fmprb_vec_init(alloc);
        TL = BL + lenbl;
        BR = TL + lentl;
        TR = BR + lenbr;
        TT = TR + lentr;

        gamma_sum_bsplit(PL, QL, BL, TL, p, q, a, m, N, L, wp);
        gamma_sum_bsplit(PR, QR, BR, TR, p, q, m, b, N, L, wp);

        /* T = BR*TL*QR + BL*TR*PL */
        _fmprb_poly_mullow(T, BR, lenbr, TL, lentl, lent, wp);
        _fmprb_vec_scalar_mul(T, T, lent, QR, wp);
        _fmprb_poly_mullow(TT, BL, lenbl, TR, lentr, lent, wp);
        _fmprb_vec_scalar_mul(TT, TT, lent, PL, wp);
        _fmprb_vec_add(T, T, TT, lent, wp);

        /* P = PL*PR, Q = QL*QR, B = BL*BR */
        fmprb_mul(P, PL, PR, wp);
        fmprb_mul(Q, QL, QR, wp);
        _fmprb_poly_mullow(B, BL, lenbl, BR, lenbr, lenb, wp);

        fmprb_clear(PL);
        fmprb_clear(QL);
        fmprb_clear(PR);
        fmprb_clear(QR);

        _fmprb_vec_clear(BL, alloc);
    }
}

static void
evaluate_series(fmprb_struct * S, const fmpq_t a, long r, long N, long len, long wp, long bsplit_wp)
{
    fmprb_t P, Q;
    fmprb_struct *B, *T;

    fmprb_init(P);
    fmprb_init(Q);

    B = _fmprb_vec_init(len);
    T = _fmprb_vec_init(len);

    gamma_sum_bsplit(P, Q, B, T, fmpq_numref(a), fmpq_denref(a), 0, r + 1, N, len, bsplit_wp);
    _fmprb_vec_scalar_mul_fmpz(T, T, len, fmpq_denref(a), bsplit_wp);

    /* S = T / (B * Q) */
    _fmprb_vec_scalar_mul(B, B, len, Q, bsplit_wp);
    _fmprb_poly_inv_series(S, B, len, bsplit_wp);
    _fmprb_poly_mullow(B, T, len, S, len, len, bsplit_wp);
    _fmprb_vec_set(S, B, len);

    fmprb_clear(P);
    fmprb_clear(Q);
    _fmprb_vec_clear(B, len);
    _fmprb_vec_clear(T, len);
}


void fmpr_gamma_ui_lbound(fmpr_t x, ulong n, long prec);

static void
add_error_series(fmprb_struct * S, long R, long N, long len)
{
    long i;
    fmpr_t t, u;

    fmpr_init(t);
    fmpr_init(u);

    fmpr_set_ui(u, N);
    fmpr_pow_sloppy_ui(t, u, R, FMPRB_RAD_PREC, FMPR_RND_UP);

    fmpr_gamma_ui_lbound(u, R, FMPRB_RAD_PREC);
    fmpr_div(t, t, u, FMPRB_RAD_PREC, FMPR_RND_UP);

    fmprb_add_error_fmpr(S + 0, t);

    for (i = 1; i < len; i++)
    {
        fmpr_div_ui(t, t, R, FMPRB_RAD_PREC, FMPR_RND_UP);
        fmprb_add_error_fmpr(S + i, t);
    }

    fmpr_clear(t);
    fmpr_clear(u);
}

static void
add_error_integral(fmprb_struct * v, long n, long num)
{
    fmpr_t t, u;
    long j;

    fmpr_init(t);
    fmpr_init(u);

    /* t = 2 exp(-n) */
    fmpr_set_si(t, -n);
    fmpr_exp(t, t, FMPRB_RAD_PREC, FMPR_RND_UP);
    fmpr_mul_2exp_si(t, t, 1);
    fmprb_add_error_fmpr(v + 0, t);

    /* u = log n */
    fmpr_set_ui(u, n);
    fmpr_log(u, u, FMPRB_RAD_PREC, FMPR_RND_UP);

    for (j = 1; j < num; j++)
    {
        /* powers of log */
        fmpr_mul(t, t, u, FMPRB_RAD_PREC, FMPR_RND_UP);
        /* factorials */
        fmpr_div_ui(t, t, j, FMPRB_RAD_PREC, FMPR_RND_UP);
        fmprb_add_error_fmpr(v + j, t);
    }

    fmpr_clear(t);
    fmpr_clear(u);
}

void
_fmprb_poly_ui_pow_series(fmprb_struct * res, ulong c, long len, long wp)
{
    long i;

    for (i = 0; i < len; i++)
    {
        if (i == 0)
        {
            fmprb_one(res + i);
        }
        else if (i == 1)
        {
            fmprb_log_ui(res + i, c, wp);
        }
        else
        {
            fmprb_mul(res + i, res + i - 1, res + 1, wp);
            fmprb_div_ui(res + i, res + i, i, wp);
        }
    }
}

void
gamma_series_fmpq_hypgeom(fmprb_struct * res, const fmpq_t a, long len, long prec)
{
    long N, R, wp, bsplit_wp;
    fmprb_struct *S, *logs;
    fmprb_t t;

    wp = prec * 1.01 + 10;
    N = wp * 0.693147180559945;

    /* round N to have many trailing zeros, speeding up arithmetic */
    if (N > 256)
    {
        int b = FLINT_BIT_COUNT(N);
        N = ((N >> (b-4)) + 1) << (b-4);
    }

    R = 3.59112147666862 * N;
    bsplit_wp = wp + N / 0.693147180559945;

    S = _fmprb_vec_init(len);
    logs = _fmprb_vec_init(len);
    fmprb_init(t);

    evaluate_series(S, a, R, N, len, wp, bsplit_wp);

    /* error from truncating the series */
    add_error_series(S, R, N, len);

    /* multiply by series for N^(a+x) */
    _fmprb_poly_ui_pow_series(logs, N, len, wp);
    _fmprb_poly_mullow(res, logs, len, S, len, len, wp);
    fmprb_set_ui(t, N);
    fmprb_pow_fmpq(t, t, a, wp);
    _fmprb_vec_scalar_mul(res, res, len, t, wp);

    /* error from truncating the integral */
    add_error_integral(res, N, len);

    _fmprb_vec_clear(S, len);
    _fmprb_vec_clear(logs, len);
    fmprb_clear(t);
}

