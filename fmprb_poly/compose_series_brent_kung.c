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
#include "fmprb_mat.h"

void
_fmprb_poly_compose_series_brent_kung(fmprb_struct * res,
    const fmprb_struct * poly1, long len1,
    const fmprb_struct * poly2, long len2, long n, long prec)
{
    fmprb_mat_t A, B, C;
    fmprb_struct *t, *h;
    long i, m;

    if (n == 1)
    {
        fmprb_set(res, poly1);
        return;
    }

    m = n_sqrt(n) + 1;

    fmprb_mat_init(A, m, n);
    fmprb_mat_init(B, m, m);
    fmprb_mat_init(C, m, n);

    h = _fmprb_vec_init(n);
    t = _fmprb_vec_init(n);

    /* Set rows of B to the segments of poly1 */
    for (i = 0; i < len1 / m; i++)
        _fmprb_vec_set(B->rows[i], poly1 + i*m, m);
    _fmprb_vec_set(B->rows[i], poly1 + i*m, len1 % m);

    /* Set rows of A to powers of poly2 */
    fmprb_set_ui(A->rows[0] + 0, 1UL);
    _fmprb_vec_set(A->rows[1], poly2, len2);
    for (i = 2; i < m; i++)
        _fmprb_poly_mullow(A->rows[i], A->rows[i-1], n, poly2, len2, n, prec);

    fmprb_mat_mul(C, B, A, prec);

    /* Evaluate block composition using the Horner scheme */
    _fmprb_vec_set(res, C->rows[m - 1], n);
    _fmprb_poly_mullow(h, A->rows[m - 1], n, poly2, len2, n, prec);

    for (i = m - 2; i >= 0; i--)
    {
        _fmprb_poly_mullow(t, res, n, h, n, n, prec);
        _fmprb_poly_add(res, t, n, C->rows[i], n, prec);
    }

    _fmprb_vec_clear(h, n);
    _fmprb_vec_clear(t, n);

    fmprb_mat_clear(A);
    fmprb_mat_clear(B);
    fmprb_mat_clear(C);
}

void
fmprb_poly_compose_series_brent_kung(fmprb_poly_t res,
                    const fmprb_poly_t poly1,
                    const fmprb_poly_t poly2, long n, long prec)
{
    long len1 = poly1->length;
    long len2 = poly2->length;
    long lenr;

    if (len2 != 0 && !fmprb_is_zero(poly2->coeffs))
    {
        printf("exception: compose_series: inner "
                "polynomial must have zero constant term\n");
        abort();
    }

    if (len1 == 0 || n == 0)
    {
        fmprb_poly_zero(res);
        return;
    }

    if (len2 == 0 || len1 == 1)
    {
        fmprb_poly_set_fmprb(res, poly1->coeffs);
        return;
    }

    lenr = FLINT_MIN((len1 - 1) * (len2 - 1) + 1, n);
    len1 = FLINT_MIN(len1, lenr);
    len2 = FLINT_MIN(len2, lenr);

    if ((res != poly1) && (res != poly2))
    {
        fmprb_poly_fit_length(res, lenr);
        _fmprb_poly_compose_series_brent_kung(res->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, prec);
        _fmprb_poly_set_length(res, lenr);
        _fmprb_poly_normalise(res);
    }
    else
    {
        fmprb_poly_t t;
        fmprb_poly_init2(t, lenr);
        _fmprb_poly_compose_series_brent_kung(t->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, prec);
        _fmprb_poly_set_length(t, lenr);
        _fmprb_poly_normalise(t);
        fmprb_poly_swap(res, t);
        fmprb_poly_clear(t);
    }
}
