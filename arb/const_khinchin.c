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
#include "arb.h"

void
arb_const_khinchin_eval_param(arb_t s, ulong N, ulong M, long prec)
{
    arb_t t, u, h;
    arb_ptr pows;
    long k, n;

    arb_init(t);
    arb_init(u);
    arb_init(h);

    if (N < 2) abort();
    /* if (M < 0) abort(); */

    pows = _arb_vec_init(N - 2);

    /* sum of logarithms */
    arb_zero(s);
    for (k = 2; k < N; k++)
    {
        arb_set_ui(t, k - 1);
        arb_div_ui(t, t, k, prec);
        arb_log(t, t, prec);

        arb_set_ui(u, k + 1);
        arb_div_ui(u, u, k, prec);
        arb_log(u, u, prec);

        arb_mul(t, t, u, prec);
        arb_sub(s, s, t, prec);
    }

    /* alternating harmonic numbers */
    arb_one(h);

    /* powers */
    for (k = 0; k < N - 2; k++)
        arb_one(pows + k);

    /* sum of zetas */
    for (n = 1; n <= M; n++)
    {
        /* zeta(2n,N) / n */
        arb_zeta_ui(t, 2 * n, prec);
        arb_sub_ui(t, t, 1, prec);
        for (k = 0; k < N - 2; k++)
        {
            arb_div_ui(pows + k, pows + k, (k + 2) * (k + 2), prec);
            arb_sub(t, t, pows + k, prec);
        }
        arb_div_ui(t, t, n, prec);

        arb_mul(t, t, h, prec);
        arb_add(s, s, t, prec);

        /* forward h by two */
        arb_set_ui(u, 2 * n);
        arb_mul_ui(u, u, 2 * n + 1, prec);
        arb_inv(u, u, prec);
        arb_sub(h, h, u, prec);
    }

    /* error bound 1/N^(2M) */
    arb_set_ui(t, N);
    arb_pow_ui(t, t, 2 * M, MAG_BITS);
    arb_inv(t, t, MAG_BITS);
    arb_add_error(s, t);

    arb_log_ui(t, 2, prec);
    arb_div(s, s, t, prec);
    arb_exp(s, s, prec);

    _arb_vec_clear(pows, N - 2);

    arb_clear(t);
    arb_clear(u);
    arb_clear(h);
}

void
arb_const_khinchin_eval(arb_t K, long prec)
{
    ulong N, M;

    prec += 10 + 2 * FLINT_BIT_COUNT(prec);

    /* heuristic */
    N = pow(prec, 0.35);

    M = (prec * 0.69314718055994530942) / (2 * log(N));

    arb_const_khinchin_eval_param(K, N, M, prec);
}

ARB_DEF_CACHED_CONSTANT(arb_const_khinchin, arb_const_khinchin_eval)
