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

    Copyright (C) 2012, 2013 Fredrik Johansson

******************************************************************************/

#include "fmprb_poly.h"

void
_fmprb_poly_log_series(fmprb_ptr res, fmprb_srcptr f, long flen, long n, long prec)
{
    flen = FLINT_MIN(flen, n);

    if (flen == 1)
    {
        fmprb_log(res, f, prec);
        _fmprb_vec_zero(res + 1, n - 1);
    }
    else if (n == 2)
    {
        fmprb_div(res + 1, f + 1, f + 0, prec);  /* safe since hlen >= 2 */
        fmprb_log(res, f, prec);
    }
    else if (_fmprb_vec_is_zero(f + 1, flen - 2))  /* f = a + bx^d */
    {
        long i, j, d = flen - 1;

        for (i = 1, j = d; j < n; j += d, i++)
        {
            if (i == 1)
                fmprb_div(res + j, f + d, f + 0, prec);
            else
                fmprb_mul(res + j, res + j - d, res + d, prec);
            _fmprb_vec_zero(res + j - d + 1, flen - 2);
        }
        _fmprb_vec_zero(res + j - d + 1, n - (j - d + 1));

        for (i = 2, j = 2 * d; j < n; j += d, i++)
            fmprb_div_si(res + j, res + j, i % 2 ? i : -i, prec);

        fmprb_log(res, f, prec); /* done last to allow aliasing */
    }
    else
    {
        fmprb_ptr f_diff, f_inv;
        fmprb_t a;
        long alloc;

        alloc = n + flen - 1;
        f_inv = _fmprb_vec_init(alloc);
        f_diff = f_inv + n;

        fmprb_init(a);
        fmprb_log(a, f, prec);

        _fmprb_poly_derivative(f_diff, f, flen, prec);
        _fmprb_poly_inv_series(f_inv, f, flen, n, prec);
        _fmprb_poly_mullow(res, f_inv, n - 1, f_diff, flen - 1, n - 1, prec);
        _fmprb_poly_integral(res, res, n, prec);
        fmprb_swap(res, a);

        fmprb_clear(a);
        _fmprb_vec_clear(f_inv, alloc);
    }
}

void
fmprb_poly_log_series(fmprb_poly_t res, const fmprb_poly_t f, long n, long prec)
{
    if (n == 0)
    {
        fmprb_poly_zero(res);
        return;
    }

    fmprb_poly_fit_length(res, n);
    if (f->length == 0)
        _fmprb_vec_indeterminate(res->coeffs, n);
    else
        _fmprb_poly_log_series(res->coeffs, f->coeffs, f->length, n, prec);
    _fmprb_poly_set_length(res, n);
    _fmprb_poly_normalise(res);
}

