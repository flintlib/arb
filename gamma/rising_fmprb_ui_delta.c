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
#include "arith.h"

void
rising_difference_polynomial(fmpz * s, fmpz * c, ulong m)
{
    long i, j, k;
    fmpz_t t;

    fmpz_init(t);

    arith_stirling_number_1u_vec(s, m, m + 1);

    for (k = 0; k < m; k++)
    {
        for (i = 0; i <= m - k - 1; i++)
        {
            fmpz_zero(c + m * k + i);

            for (j = i + 1; j + k <= m; j++)
            {
                if (j == i + 1)
                {
                    fmpz_bin_uiui(t, 1+i+k, k);
                    fmpz_mul_ui(t, t, m * (i+1));
                }
                else
                {
                    fmpz_mul_ui(t, t, m * (k + j));
                    fmpz_divexact_ui(t, t, j - i);
                }

                fmpz_addmul(c + m * k + i, s + j + k, t);
            }
        }
    }

    fmpz_clear(t);
}

void
gamma_rising_fmprb_ui_delta(fmprb_t y, const fmprb_t x, ulong n, ulong m, long prec)
{
    fmprb_struct * xs;
    fmprb_t t, u, v;
    ulong i, k, rem;
    fmpz_t c, h;
    fmpz *s, *d;
    long wp;

    if (n == 0)
    {
        fmprb_one(y);
        return;
    }

    if (n == 1)
    {
        fmprb_set_round(y, x, prec);
        return;
    }

    wp = FMPR_PREC_ADD(prec, FLINT_BIT_COUNT(n));

    fmprb_init(t);
    fmprb_init(u);
    fmprb_init(v);
    fmpz_init(c);
    fmpz_init(h);

    if (m == 0)
    {
        ulong m1, m2;
        m1 = 0.2 * pow(wp, 0.4);
        m2 = n_sqrt(n);
        m = FLINT_MIN(m1, m2);
    }

    m = FLINT_MIN(m, n);
    m = FLINT_MAX(m, 1);

    xs = _fmprb_vec_init(m + 1);
    d = _fmpz_vec_init(m * m);
    s = _fmpz_vec_init(m + 1);

    for (i = 0; i <= m; i++)
    {
        if (i == 0)
            fmprb_one(xs + i);
        else if (i == 1)
            fmprb_set(xs + i, x);
        else if (i % 2 == 0)
            fmprb_mul(xs + i, xs + i / 2, xs + i / 2, wp);
        else
            fmprb_mul(xs + i, xs + i - 1, x, wp);
    }

    rising_difference_polynomial(s, d, m);

    /* tail */
    rem = m;
    while (rem + m <= n)
        rem += m;
    fmprb_one(y);
    for (k = rem; k < n; k++)
    {
        fmprb_add_ui(t, xs + 1, k, wp);
        fmprb_mul(y, y, t, wp);
    }

    /* initial rising factorial */
    fmprb_zero(t);
    for (i = 1; i <= m; i++)
        fmprb_addmul_fmpz(t, xs + i, s + i, wp);

    fmprb_mul(y, y, t, wp);

    /* the leading coefficient is always the same */
    fmprb_mul_fmpz(xs + m - 1, xs + m - 1, d + m - 1 + 0, wp);

    for (k = 0; k + 2 * m <= n; k += m)
    {
        for (i = 0; i < m - 1; i++)
        {
            fmpz_set_ui(h, k);
            _fmpz_poly_evaluate_horner_fmpz(c, d + i * m, m - i, h);

            if (i == 0)
                fmprb_add_fmpz(t, t, c, wp);
            else
                fmprb_addmul_fmpz(t, xs + i, c, wp);
        }

        fmprb_add(t, t, xs + m - 1, wp);
        fmprb_mul(y, y, t, wp);
    }

    fmprb_set_round(y, y, prec);

    fmprb_clear(t);
    fmprb_clear(u);
    fmprb_clear(v);
    _fmprb_vec_clear(xs, m + 1);
    _fmpz_vec_clear(d, m * m);
    _fmpz_vec_clear(s, m + 1);
    fmpz_clear(c);
    fmpz_clear(h);
}

