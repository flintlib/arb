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

#include "arb_poly.h"

/* Bound based on binomial theorem */
long
_arb_poly_swinnerton_dyer_ui_prec(ulong n)
{
    long i;
    double u, N;

    N = 1UL << n;

    /* u = (sum of square roots)^(2^n) */
    u = 0;
    for (i = 0; i < n; i++)
        u += sqrt(n_nth_prime(1 + i));
    u = N * log(u) * 1.44269504088897;

    /* Central binomial coefficient C(N,N/2) < 2^N / sqrt(3*N/2) */
    u += N - 0.5*(n-1) - 0.792481250360578; /* log(sqrt(3)) */

    /* experimental heuristic: the bound is 2x too large */
    return u * 0.5 + 15;
}

void
_arb_poly_swinnerton_dyer_ui(arb_ptr T, ulong n, long trunc, long prec)
{
    arb_ptr square_roots, tmp1, tmp2, tmp3;
    arb_t one;
    long i, j, k, N;

    if (n == 0)
    {
        arb_zero(T);
        arb_one(T + 1);
        return;
    }

    if (prec == 0)
        prec = _arb_poly_swinnerton_dyer_ui_prec(n);

    N = 1L << n;
    trunc = FLINT_MIN(trunc, N + 1);

    arb_init(one);
    arb_one(one);

    square_roots = _arb_vec_init(n);
    tmp1 = flint_malloc((N/2 + 1) * sizeof(arb_struct));
    tmp2 = flint_malloc((N/2 + 1) * sizeof(arb_struct));
    tmp3 = _arb_vec_init(N);

    for (i = 0; i < n; i++)
        arb_sqrt_ui(square_roots + i, n_nth_prime(i + 1), prec);

    /* Build linear factors */
    for (i = 0; i < N; i++)
    {
        arb_zero(T + i);

        for (j = 0; j < n; j++)
        {
            if ((i >> j) & 1)
                arb_add(T + i, T + i, square_roots + j, prec);
            else
                arb_sub(T + i, T + i, square_roots + j, prec);
        }
    }

    /* For each level... */
    for (i = 0; i < n; i++)
    {
        long stride = 1UL << i;

        for (j = 0; j < N; j += 2*stride)
        {
            for (k = 0; k < stride; k++)
            {
                tmp1[k] = T[j + k];
                tmp2[k] = T[j + stride + k];
            }

            tmp1[stride] = *one;
            tmp2[stride] = *one;

            _arb_poly_mullow(tmp3, tmp1, stride + 1, tmp2, stride + 1,
                FLINT_MIN(2 * stride, trunc), prec);
            _arb_vec_set(T + j, tmp3, FLINT_MIN(2 * stride, trunc));
        }
    }

    arb_one(T + N);
    _arb_vec_clear(square_roots, n);
    flint_free(tmp1);
    flint_free(tmp2);
    _arb_vec_clear(tmp3, 1UL << n);
    arb_clear(one);
}

void
arb_poly_swinnerton_dyer_ui(arb_poly_t poly, ulong n, long prec)
{
    long N = 1L << n;

    arb_poly_fit_length(poly, N + 1);
    _arb_poly_swinnerton_dyer_ui(poly->coeffs, n, N + 1, prec);
    _arb_poly_set_length(poly, N + 1);
    _arb_poly_normalise(poly);
}

