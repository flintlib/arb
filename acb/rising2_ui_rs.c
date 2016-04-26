/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "flint/fmpz_poly.h"
#include "acb.h"

void _gamma_rf_bsplit(fmpz * A, ulong a, ulong b);

void
acb_rising2_ui_rs(acb_t u, acb_t v,
    const acb_t x, ulong n, ulong m, slong prec)
{
    if (n == 0)
    {
        acb_zero(v);
        acb_one(u);
    }
    else if (n == 1)
    {
        acb_set(u, x);
        acb_one(v);
    }
    else
    {
        slong wp;
        ulong i, j, a, b;
        acb_ptr xs;
        acb_t S, T, U, V;
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

        xs = _acb_vec_init(m + 1);
        A = _fmpz_vec_init(2 * m + 1);
        B = A + (m + 1);

        acb_init(S);
        acb_init(T);
        acb_init(U);
        acb_init(V);
        _acb_vec_set_powers(xs, x, m + 1, wp);

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

            acb_set_fmpz(S, A);

            for (j = 1; j <= b - a; j++)
                acb_addmul_fmpz(S, xs + j, A + j, wp);

            acb_set_fmpz(T, B);

            for (j = 1; j < b - a; j++)
                acb_addmul_fmpz(T, xs + j, B + j, wp);

            if (i == 0)
            {
                acb_set(U, S);
                acb_set(V, T);
            }
            else
            {
                acb_mul(V, V, S, wp);
                acb_addmul(V, U, T, wp);
                acb_mul(U, U, S, wp);
            }
        }

        acb_set(u, U);
        acb_set(v, V);

        _acb_vec_clear(xs, m + 1);
        _fmpz_vec_clear(A, 2 * m + 1);

        acb_clear(S);
        acb_clear(T);
        acb_clear(U);
        acb_clear(V);
    }
}

