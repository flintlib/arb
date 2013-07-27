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

#include "fmpcb_poly.h"

void
_fmpcb_poly_compose_series_horner(fmpcb_ptr res, fmpcb_srcptr poly1, long len1,
                            fmpcb_srcptr poly2, long len2, long n, long prec)
{
    if (n == 1)
    {
        fmpcb_set(res, poly1);
    }
    else
    {
        long i = len1 - 1;
        long lenr;

        fmpcb_ptr t = _fmpcb_vec_init(n);

        lenr = len2;
        _fmpcb_vec_scalar_mul(res, poly2, len2, poly1 + i, prec);
        i--;
        fmpcb_add(res, res, poly1 + i, prec);

        while (i > 0)
        {
            i--;
            if (lenr + len2 - 1 < n)
            {
                _fmpcb_poly_mul(t, res, lenr, poly2, len2, prec);
                lenr = lenr + len2 - 1;
            }
            else
            {
                _fmpcb_poly_mullow(t, res, lenr, poly2, len2, n, prec);
                lenr = n;
            }
            _fmpcb_poly_add(res, t, lenr, poly1 + i, 1, prec);
        }

        _fmpcb_vec_zero(res + lenr, n - lenr);
        _fmpcb_vec_clear(t, n);
    }
}

void
fmpcb_poly_compose_series_horner(fmpcb_poly_t res,
                    const fmpcb_poly_t poly1,
                    const fmpcb_poly_t poly2, long n, long prec)
{
    long len1 = poly1->length;
    long len2 = poly2->length;
    long lenr;

    if (len2 != 0 && !fmpcb_is_zero(poly2->coeffs))
    {
        printf("exception: compose_series: inner "
                "polynomial must have zero constant term\n");
        abort();
    }

    if (len1 == 0 || n == 0)
    {
        fmpcb_poly_zero(res);
        return;
    }

    if (len2 == 0 || len1 == 1)
    {
        fmpcb_poly_set_fmpcb(res, poly1->coeffs);
        return;
    }

    lenr = FLINT_MIN((len1 - 1) * (len2 - 1) + 1, n);
    len1 = FLINT_MIN(len1, lenr);
    len2 = FLINT_MIN(len2, lenr);

    if ((res != poly1) && (res != poly2))
    {
        fmpcb_poly_fit_length(res, lenr);
        _fmpcb_poly_compose_series_horner(res->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, prec);
        _fmpcb_poly_set_length(res, lenr);
        _fmpcb_poly_normalise(res);
    }
    else
    {
        fmpcb_poly_t t;
        fmpcb_poly_init2(t, lenr);
        _fmpcb_poly_compose_series_horner(t->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, prec);
        _fmpcb_poly_set_length(t, lenr);
        _fmpcb_poly_normalise(t);
        fmpcb_poly_swap(res, t);
        fmpcb_poly_clear(t);
    }
}
