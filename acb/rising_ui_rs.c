/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "flint/fmpz_poly.h"
#include "acb.h"

void rising_difference_polynomial(fmpz * s, fmpz * c, ulong m);

void
acb_rising_ui_rs(acb_t y, const acb_t x, ulong n, ulong m, slong prec)
{
    acb_ptr xs;
    acb_t t, u, v;
    ulong i, k, rem;
    fmpz_t c, h;
    fmpz *s, *d;
    slong wp;

    if (n == 0)
    {
        acb_one(y);
        return;
    }

    if (n == 1)
    {
        acb_set_round(y, x, prec);
        return;
    }

    wp = ARF_PREC_ADD(prec, FLINT_BIT_COUNT(n));

    acb_init(t);
    acb_init(u);
    acb_init(v);
    fmpz_init(c);
    fmpz_init(h);

    if (m == 0)
    {
        ulong m1, m2;
        m1 = 0.2 * pow(2.0 * wp, 0.4);
        m2 = n_sqrt(n);
        m = FLINT_MIN(m1, m2);
    }

    m = FLINT_MIN(m, n);
    m = FLINT_MAX(m, 1);

    xs = _acb_vec_init(m + 1);
    d = _fmpz_vec_init(m * m);
    s = _fmpz_vec_init(m + 1);

    _acb_vec_set_powers(xs, x, m + 1, wp);

    rising_difference_polynomial(s, d, m);

    /* tail */
    rem = m;
    while (rem + m <= n)
        rem += m;
    acb_one(y);
    for (k = rem; k < n; k++)
    {
        acb_add_ui(t, xs + 1, k, wp);
        acb_mul(y, y, t, wp);
    }

    /* initial rising factorial */
    acb_zero(t);
    for (i = 1; i <= m; i++)
        acb_addmul_fmpz(t, xs + i, s + i, wp);

    acb_mul(y, y, t, wp);

    /* the leading coefficient is always the same */
    acb_mul_fmpz(xs + m - 1, xs + m - 1, d + m - 1 + 0, wp);

    for (k = 0; k + 2 * m <= n; k += m)
    {
        for (i = 0; i < m - 1; i++)
        {
            fmpz_set_ui(h, k);
            _fmpz_poly_evaluate_horner_fmpz(c, d + i * m, m - i, h);

            if (i == 0)
                acb_add_fmpz(t, t, c, wp);
            else
                acb_addmul_fmpz(t, xs + i, c, wp);
        }

        acb_add(t, t, xs + m - 1, wp);
        acb_mul(y, y, t, wp);
    }

    acb_set_round(y, y, prec);

    acb_clear(t);
    acb_clear(u);
    acb_clear(v);
    _acb_vec_clear(xs, m + 1);
    _fmpz_vec_clear(d, m * m);
    _fmpz_vec_clear(s, m + 1);
    fmpz_clear(c);
    fmpz_clear(h);
}

