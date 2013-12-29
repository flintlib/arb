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

#include "gamma.h"

void _gamma_rf_bsplit(fmpz * A, ulong a, ulong b);

void
gamma_rising2_fmpcb_ui_rs(fmpcb_t u, fmpcb_t v,
    const fmpcb_t x, ulong n, ulong m, long prec)
{
    if (n == 0)
    {
        fmpcb_zero(v);
        fmpcb_one(u);
    }
    else if (n == 1)
    {
        fmpcb_set(u, x);
        fmpcb_one(v);
    }
    else
    {
        long wp;
        ulong i, j, a, b;
        fmpcb_ptr xs;
        fmpcb_t S, T, U, V;
        fmpz *A, *B;

        wp = FMPR_PREC_ADD(prec, FLINT_BIT_COUNT(n));

        if (m == 0)
        {
            ulong m1, m2;
            m1 = 0.6 * pow(wp, 0.4);
            m2 = n_sqrt(n);
            m = FLINT_MIN(m1, m2);
        }

        m = FLINT_MAX(m, 1);

        xs = _fmpcb_vec_init(m + 1);
        A = _fmpz_vec_init(2 * m + 1);
        B = A + (m + 1);

        fmpcb_init(S);
        fmpcb_init(T);
        fmpcb_init(U);
        fmpcb_init(V);
        _fmpcb_vec_set_powers(xs, x, m + 1, wp);

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

            fmpcb_set_fmpz(S, A);

            for (j = 1; j <= b - a; j++)
                fmpcb_addmul_fmpz(S, xs + j, A + j, wp);

            fmpcb_set_fmpz(T, B);

            for (j = 1; j < b - a; j++)
                fmpcb_addmul_fmpz(T, xs + j, B + j, wp);

            if (i == 0)
            {
                fmpcb_set(U, S);
                fmpcb_set(V, T);
            }
            else
            {
                fmpcb_mul(V, V, S, wp);
                fmpcb_addmul(V, U, T, wp);
                fmpcb_mul(U, U, S, wp);
            }
        }

        fmpcb_set(u, U);
        fmpcb_set(v, V);

        _fmpcb_vec_clear(xs, m + 1);
        _fmpz_vec_clear(A, 2 * m + 1);

        fmpcb_clear(S);
        fmpcb_clear(T);
        fmpcb_clear(U);
        fmpcb_clear(V);
    }
}

