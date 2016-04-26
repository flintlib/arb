/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "flint/fmpz_poly.h"
#include "arb.h"

void
_gamma_rf_bsplit(fmpz * A, ulong a, ulong b)
{
    ulong n = b - a;

    if (n == 0)
    {
        fmpz_one(A);
    }
    else if (n < 8)
    {
        ulong j, k;

        fmpz_set_ui(A, a);
        fmpz_one(A + 1);

        for (j = 1; j < n; j++)
        {
            fmpz_one(A + j + 1);

            for (k = j; k > 0; k--)
            {
                fmpz_mul_ui(A + k, A + k, a + j);
                fmpz_add(A + k, A + k, A + k - 1);
            }

            fmpz_mul_ui(A, A, a + j);
        }
    }
    else
    {
        ulong m = a + (b - a) / 2;
        ulong w = m - a;
        ulong v = b - m;

        fmpz *t, *A1, *A2;

        t = _fmpz_vec_init(w + v + 2);

        A1 = t;
        A2 = A1 + w + 1;

        _gamma_rf_bsplit(A1, a, m);
        _gamma_rf_bsplit(A2, m, b);

        _fmpz_poly_mul(A, A2, v + 1, A1, w + 1);

        _fmpz_vec_clear(t, w + v + 2);
    }
}

void
arb_rising2_ui_rs(arb_t u, arb_t v,
    const arb_t x, ulong n, ulong m, slong prec)
{
    if (n == 0)
    {
        arb_zero(v);
        arb_one(u);
    }
    else if (n == 1)
    {
        arb_set(u, x);
        arb_one(v);
    }
    else
    {
        slong wp;
        ulong i, j, a, b;
        arb_ptr xs;
        arb_t S, T, U, V;
        fmpz *A, *B;

        wp = ARF_PREC_ADD(prec, FLINT_BIT_COUNT(n));

        if (m == 0)
        {
            ulong m1, m2;
            m1 = 0.6 * pow(wp, 0.4);
            m2 = n_sqrt(n);
            m = FLINT_MIN(m1, m2);
        }

        m = FLINT_MAX(m, 1);

        xs = _arb_vec_init(m + 1);
        A = _fmpz_vec_init(2 * m + 1);
        B = A + (m + 1);

        arb_init(S);
        arb_init(T);
        arb_init(U);
        arb_init(V);
        _arb_vec_set_powers(xs, x, m + 1, wp);

        for (i = 0; i < n; i += m)
        {
            a = i;
            b = FLINT_MIN(n, a + m);

            if (a == 0 || b != a + m)
            {
                _gamma_rf_bsplit(A, a, b);
            }
            else
            {
                fmpz tt = m;
                _fmpz_poly_taylor_shift(A, &tt, m + 1);
            }

            _fmpz_poly_derivative(B, A, b - a + 1);

            arb_set_fmpz(S, A);

            for (j = 1; j <= b - a; j++)
                arb_addmul_fmpz(S, xs + j, A + j, wp);

            arb_set_fmpz(T, B);

            for (j = 1; j < b - a; j++)
                arb_addmul_fmpz(T, xs + j, B + j, wp);

            if (i == 0)
            {
                arb_set(U, S);
                arb_set(V, T);
            }
            else
            {
                arb_mul(V, V, S, wp);
                arb_addmul(V, U, T, wp);
                arb_mul(U, U, S, wp);
            }
        }

        arb_set(u, U);
        arb_set(v, V);

        _arb_vec_clear(xs, m + 1);
        _fmpz_vec_clear(A, 2 * m + 1);

        arb_clear(S);
        arb_clear(T);
        arb_clear(U);
        arb_clear(V);
    }
}

