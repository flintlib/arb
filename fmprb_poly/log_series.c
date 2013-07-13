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

#include "fmprb_poly.h"

void
_fmprb_poly_log_series(fmprb_struct * res, fmprb_struct * f, long n, long prec)
{
    fmprb_struct * f_diff;
    fmprb_struct * f_inv;
    fmprb_t a;

    f_diff = _fmprb_vec_init(n);
    f_inv = _fmprb_vec_init(n);
    fmprb_init(a);

    fmprb_log(a, f, prec);

    _fmprb_poly_derivative(f_diff, f, n, prec);
    fmprb_zero(f_diff + n - 1);
    _fmprb_poly_inv_series(f_inv, f, n, n, prec);
    _fmprb_poly_mullow(res, f_diff, n - 1, f_inv, n - 1, n - 1, prec);
    _fmprb_poly_integral(res, res, n, prec);
    fmprb_set(res, a);

    fmprb_clear(a);
    _fmprb_vec_clear(f_diff, n);
    _fmprb_vec_clear(f_inv, n);
}

void
fmprb_poly_log_series(fmprb_poly_t res, const fmprb_poly_t f, long n, long prec)
{
    if (f->length < 1 || n < 1)
    {
        printf("Exception: log_series: input must be nonzero\n");
        abort();
    }

    if (f->length < n || res == f)
    {
        long i;
        fmprb_poly_t t;

        fmprb_poly_init(t);
        fmprb_poly_fit_length(t, n);

        for (i = 0; i < FLINT_MIN(n, f->length); i++)
            fmprb_set(t->coeffs + i, f->coeffs + i);
        for (i = FLINT_MIN(n, f->length); i < n; i++)
            fmprb_zero(t->coeffs + i);

        _fmprb_poly_set_length(t, n);
        fmprb_poly_log_series(res, t, n, prec);

        fmprb_poly_clear(t);
        return;
    }

    fmprb_poly_fit_length(res, n);
    _fmprb_poly_log_series(res->coeffs, f->coeffs, n, prec);
    _fmprb_poly_set_length(res, n);
    _fmprb_poly_normalise(res);
}
