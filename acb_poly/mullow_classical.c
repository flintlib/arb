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

    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "acb_poly.h"

void
_acb_poly_mullow_classical(acb_ptr res,
    acb_srcptr poly1, long len1,
    acb_srcptr poly2, long len2, long n, long prec)
{
    len1 = FLINT_MIN(len1, n);
    len2 = FLINT_MIN(len2, n);

    if (n == 1)
    {
        acb_mul(res, poly1, poly2, prec);
    }
    else if (poly1 == poly2 && len1 == len2)
    {
        long i;

        _acb_vec_scalar_mul(res, poly1, FLINT_MIN(len1, n), poly1, prec);
        _acb_vec_scalar_mul(res + len1, poly1 + 1, n - len1, poly1 + len1 - 1, prec);

        for (i = 1; i < len1 - 1; i++)
            _acb_vec_scalar_addmul(res + i + 1, poly1 + 1,
                FLINT_MIN(i - 1, n - (i + 1)), poly1 + i, prec);

        for (i = 1; i < FLINT_MIN(2 * len1 - 2, n); i++)
            acb_mul_2exp_si(res + i, res + i, 1);

        for (i = 1; i < FLINT_MIN(len1 - 1, (n + 1) / 2); i++)
            acb_addmul(res + 2 * i, poly1 + i, poly1 + i, prec);
    }
    else
    {
        long i;

        _acb_vec_scalar_mul(res, poly1, FLINT_MIN(len1, n), poly2, prec);

        if (n > len1)
            _acb_vec_scalar_mul(res + len1, poly2 + 1, n - len1,
                                      poly1 + len1 - 1, prec);

        for (i = 0; i < FLINT_MIN(len1, n) - 1; i++)
            _acb_vec_scalar_addmul(res + i + 1, poly2 + 1,
                                         FLINT_MIN(len2, n - i) - 1,
                                         poly1 + i, prec);
    }
}

void
acb_poly_mullow_classical(acb_poly_t res, const acb_poly_t poly1,
                                            const acb_poly_t poly2,
                                                long n, long prec)
{
    long len_out;

    if (poly1->length == 0 || poly2->length == 0 || n == 0)
    {
        acb_poly_zero(res);
        return;
    }

    len_out = poly1->length + poly2->length - 1;
    if (n > len_out)
        n = len_out;

    if (res == poly1 || res == poly2)
    {
        acb_poly_t t;
        acb_poly_init2(t, n);
        _acb_poly_mullow_classical(t->coeffs, poly1->coeffs, poly1->length,
                                    poly2->coeffs, poly2->length, n, prec);
        acb_poly_swap(res, t);
        acb_poly_clear(t);
    }
    else
    {
        acb_poly_fit_length(res, n);
        _acb_poly_mullow_classical(res->coeffs, poly1->coeffs, poly1->length,
                                    poly2->coeffs, poly2->length, n, prec);
    }

    _acb_poly_set_length(res, n);
    _acb_poly_normalise(res);
}

