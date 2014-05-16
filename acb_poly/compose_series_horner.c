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

#include "acb_poly.h"

void
_acb_poly_compose_series_horner(acb_ptr res, acb_srcptr poly1, long len1,
                            acb_srcptr poly2, long len2, long n, long prec)
{
    if (n == 1)
    {
        acb_set(res, poly1);
    }
    else
    {
        long i = len1 - 1;
        long lenr;

        acb_ptr t = _acb_vec_init(n);

        lenr = len2;
        _acb_vec_scalar_mul(res, poly2, len2, poly1 + i, prec);
        i--;
        acb_add(res, res, poly1 + i, prec);

        while (i > 0)
        {
            i--;
            if (lenr + len2 - 1 < n)
            {
                _acb_poly_mul(t, res, lenr, poly2, len2, prec);
                lenr = lenr + len2 - 1;
            }
            else
            {
                _acb_poly_mullow(t, res, lenr, poly2, len2, n, prec);
                lenr = n;
            }
            _acb_poly_add(res, t, lenr, poly1 + i, 1, prec);
        }

        _acb_vec_zero(res + lenr, n - lenr);
        _acb_vec_clear(t, n);
    }
}

void
acb_poly_compose_series_horner(acb_poly_t res,
                    const acb_poly_t poly1,
                    const acb_poly_t poly2, long n, long prec)
{
    long len1 = poly1->length;
    long len2 = poly2->length;
    long lenr;

    if (len2 != 0 && !acb_is_zero(poly2->coeffs))
    {
        printf("exception: compose_series: inner "
                "polynomial must have zero constant term\n");
        abort();
    }

    if (len1 == 0 || n == 0)
    {
        acb_poly_zero(res);
        return;
    }

    if (len2 == 0 || len1 == 1)
    {
        acb_poly_set_acb(res, poly1->coeffs);
        return;
    }

    lenr = FLINT_MIN((len1 - 1) * (len2 - 1) + 1, n);
    len1 = FLINT_MIN(len1, lenr);
    len2 = FLINT_MIN(len2, lenr);

    if ((res != poly1) && (res != poly2))
    {
        acb_poly_fit_length(res, lenr);
        _acb_poly_compose_series_horner(res->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, prec);
        _acb_poly_set_length(res, lenr);
        _acb_poly_normalise(res);
    }
    else
    {
        acb_poly_t t;
        acb_poly_init2(t, lenr);
        _acb_poly_compose_series_horner(t->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, prec);
        _acb_poly_set_length(t, lenr);
        _acb_poly_normalise(t);
        acb_poly_swap(res, t);
        acb_poly_clear(t);
    }
}
