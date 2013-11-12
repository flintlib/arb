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

#include <math.h>
#include "zeta.h"

void
fmprb_const_khinchin_eval_param(fmprb_t s, ulong N, ulong M, long prec)
{
    fmprb_t t, u, h;
    fmprb_ptr pows;
    long k, n;

    fmprb_init(t);
    fmprb_init(u);
    fmprb_init(h);

    if (N < 2) abort();
    if (M < 0) abort();

    pows = _fmprb_vec_init(N - 2);

    /* sum of logarithms */
    fmprb_zero(s);
    for (k = 2; k < N; k++)
    {
        fmprb_set_ui(t, k - 1);
        fmprb_div_ui(t, t, k, prec);
        fmprb_log(t, t, prec);

        fmprb_set_ui(u, k + 1);
        fmprb_div_ui(u, u, k, prec);
        fmprb_log(u, u, prec);

        fmprb_mul(t, t, u, prec);
        fmprb_sub(s, s, t, prec);
    }

    /* alternating harmonic numbers */
    fmprb_one(h);

    /* powers */
    for (k = 0; k < N - 2; k++)
        fmprb_one(pows + k);

    /* sum of zetas */
    for (n = 1; n <= M; n++)
    {
        /* zeta(2n,N) / n */
        zeta_ui(t, 2 * n, prec);
        fmprb_sub_ui(t, t, 1, prec);
        for (k = 0; k < N - 2; k++)
        {
            fmprb_div_ui(pows + k, pows + k, (k + 2) * (k + 2), prec);
            fmprb_sub(t, t, pows + k, prec);
        }
        fmprb_div_ui(t, t, n, prec);

        fmprb_mul(t, t, h, prec);
        fmprb_add(s, s, t, prec);

        /* forward h by two */
        fmprb_set_ui(u, 2 * n);
        fmprb_mul_ui(u, u, 2 * n + 1, prec);
        fmprb_inv(u, u, prec);
        fmprb_sub(h, h, u, prec);
    }

    /* error bound 1/N^(2M) */
    fmprb_set_ui(t, N);
    fmprb_pow_ui(t, t, 2 * M, FMPRB_RAD_PREC);
    fmprb_inv(t, t, FMPRB_RAD_PREC);
    fmprb_add_error(s, t);

    fmprb_log_ui(t, 2, prec);
    fmprb_div(s, s, t, prec);
    fmprb_exp(s, s, prec);

    _fmprb_vec_clear(pows, N - 2);

    fmprb_clear(t);
    fmprb_clear(u);
    fmprb_clear(h);
}

void
fmprb_const_khinchin_eval(fmprb_t K, long prec)
{
    ulong N, M;

    prec += 10 + 2 * FLINT_BIT_COUNT(prec);

    /* heuristic */
    N = pow(prec, 0.35);

    M = (prec * 0.69314718055994530942) / (2 * log(N));

    fmprb_const_khinchin_eval_param(K, N, M, prec);
}

DEF_CACHED_CONSTANT(fmprb_const_khinchin, fmprb_const_khinchin_eval)

