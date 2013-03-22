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

void
gamma_fmpq_outward(fmprb_t y, const fmpq_t x, long prec)
{
    fmpq_t a;
    fmpz_t n;
    fmpz p, q;
    long m;
    fmprb_t t, u;

    fmpq_init(a);
    fmpz_init(n);
    fmprb_init(t);
    fmprb_init(u);

    /* write x = a + n with 0 < a <= 1 */
    if (fmpz_is_one(fmpq_denref(x)))
    {
        fmpq_one(a);
        fmpz_sub_ui(n, fmpq_numref(x), 1);
    }
    else
    {
        fmpz_fdiv_qr(n, fmpq_numref(a), fmpq_numref(x), fmpq_denref(x));
        fmpz_set(fmpq_denref(a), fmpq_denref(x));
    }

    if (!fmpz_fits_si(n))
    {
        printf("gamma: too large fmpq to reduce to 0!\n");
        abort();
    }

    m = fmpz_get_si(n);

    /* evaluate gamma(a) */
    p = *fmpq_numref(a);
    q = *fmpq_denref(a);

    if (q == 1 || q == 2 || q == 3 || q == 4 || q == 6)
        gamma_small_frac(t, p, q, prec);
    else
        gamma_series_fmpq_hypgeom(t, a, 1, prec);

    /* argument reduction */
    if (m >= 0)
    {
        gamma_rising_fmprb_fmpq_ui_bsplit(u, a, m, prec);
        fmprb_mul(y, t, u, prec);
    }
    else
    {
        gamma_rising_fmprb_fmpq_ui_bsplit(u, x, -m, prec);
        fmprb_div(y, t, u, prec);
    }

    fmpq_clear(a);
    fmpz_clear(n);
    fmprb_clear(t);
    fmprb_clear(u);
}

