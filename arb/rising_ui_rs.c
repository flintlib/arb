/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "flint/arith.h"
#include "arb.h"

void
rising_difference_polynomial(fmpz * s, fmpz * c, ulong m)
{
    slong i, j, v;

    fmpz_t t;
    fmpz_init(t);

    arith_stirling_number_1u_vec(s, m, m + 1);

    /* Compute the first row */
    for (i = 0; i < m; i++)
    {
        fmpz_set_ui(t, m);
        fmpz_mul_ui(t, t, (i + 1));
        fmpz_mul(c + i, s + (i + 1), t);

        for (j = i + 2; j < m + 1; j++)
        {
            fmpz_mul_ui(t, t, m * j);
            fmpz_divexact_ui(t, t, j - i);
            fmpz_addmul(c + i, s + j, t);
        }
    }

    /* Extend using recurrence and symmetry */
    for (v = 1; v < m; v++)
    {
        for (i = v; i < m - v; i++)
        {
            fmpz_mul_ui(t, c + (v - 1) * m + (i + 1), i + 1);
            fmpz_divexact_ui(c + v * m + i, t, v);
        }

        for (i = 0; i < v; i++)
            fmpz_set(c + v * m + i, c + i * m + v);
    }

    fmpz_clear(t);
}

void
arb_rising_ui_rs(arb_t y, const arb_t x, ulong n, ulong m, slong prec)
{
    arb_ptr xs;
    arb_t t, u, v;
    ulong i, k, rem;
    fmpz_t c, h;
    fmpz *s, *d;
    slong wp;

    if (n == 0)
    {
        arb_one(y);
        return;
    }

    if (n == 1)
    {
        arb_set_round(y, x, prec);
        return;
    }

    wp = ARF_PREC_ADD(prec, FLINT_BIT_COUNT(n));

    arb_init(t);
    arb_init(u);
    arb_init(v);
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

    xs = _arb_vec_init(m + 1);
    d = _fmpz_vec_init(m * m);
    s = _fmpz_vec_init(m + 1);

    _arb_vec_set_powers(xs, x, m + 1, wp);

    rising_difference_polynomial(s, d, m);

    /* tail */
    rem = m;
    while (rem + m <= n)
        rem += m;
    arb_one(y);
    for (k = rem; k < n; k++)
    {
        arb_add_ui(t, xs + 1, k, wp);
        arb_mul(y, y, t, wp);
    }

    /* initial rising factorial */
    arb_zero(t);
    for (i = 1; i <= m; i++)
        arb_addmul_fmpz(t, xs + i, s + i, wp);

    arb_mul(y, y, t, wp);

    /* the leading coefficient is always the same */
    arb_mul_fmpz(xs + m - 1, xs + m - 1, d + m - 1 + 0, wp);

    for (k = 0; k + 2 * m <= n; k += m)
    {
        for (i = 0; i < m - 1; i++)
        {
            fmpz_set_ui(h, k);
            _fmpz_poly_evaluate_horner_fmpz(c, d + i * m, m - i, h);

            if (i == 0)
                arb_add_fmpz(t, t, c, wp);
            else
                arb_addmul_fmpz(t, xs + i, c, wp);
        }

        arb_add(t, t, xs + m - 1, wp);
        arb_mul(y, y, t, wp);
    }

    arb_set_round(y, y, prec);

    arb_clear(t);
    arb_clear(u);
    arb_clear(v);
    _arb_vec_clear(xs, m + 1);
    _fmpz_vec_clear(d, m * m);
    _fmpz_vec_clear(s, m + 1);
    fmpz_clear(c);
    fmpz_clear(h);
}

