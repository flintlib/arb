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
_fmpcb_poly_mullow_transpose(fmpcb_struct * res,
    const fmpcb_struct * poly1, long len1,
    const fmpcb_struct * poly2, long len2, long n, long prec)
{
    fmprb_struct *a, *b, *c, *d, *e, *f, *w;
    fmprb_struct *t, *u, *v;
    long i;

    len1 = FLINT_MIN(len1, n);
    len2 = FLINT_MIN(len2, n);

    w = flint_malloc(sizeof(fmprb_struct) * (2 * (len1 + len2 + n)));
    a = w;
    b = a + len1;
    c = b + len1;
    d = c + len2;
    e = d + len2;
    f = e + n;

    t = _fmprb_vec_init(n);
    u = _fmprb_vec_init(n);
    v = _fmprb_vec_init(n);

    for (i = 0; i < len1; i++)
    {
        a[i] = *fmpcb_realref(poly1 + i);
        b[i] = *fmpcb_imagref(poly1 + i);
    }

    for (i = 0; i < len2; i++)
    {
        c[i] = *fmpcb_realref(poly2 + i);
        d[i] = *fmpcb_imagref(poly2 + i);
    }

    for (i = 0; i < n; i++)
    {
        e[i] = *fmpcb_realref(res + i);
        f[i] = *fmpcb_imagref(res + i);
    }

    _fmprb_vec_add(t, a, b, len1, prec);
    _fmprb_vec_add(u, c, d, len2, prec);

    _fmprb_poly_mullow(v, t, len1, u, len2, n, prec);
    _fmprb_poly_mullow(t, a, len1, c, len2, n, prec);
    _fmprb_poly_mullow(u, b, len1, d, len2, n, prec);

    _fmprb_vec_sub(e, t, u, n, prec);
    _fmprb_vec_sub(f, v, t, n, prec);
    _fmprb_vec_sub(f, f, u, n, prec);

    for (i = 0; i < n; i++)
    {
        *fmpcb_realref(res + i) = e[i];
        *fmpcb_imagref(res + i) = f[i];
    }

    _fmprb_vec_clear(t, n);
    _fmprb_vec_clear(u, n);
    _fmprb_vec_clear(v, n);

    flint_free(w);
}

void
fmpcb_poly_mullow_transpose(fmpcb_poly_t res, const fmpcb_poly_t poly1,
                                            const fmpcb_poly_t poly2,
                                                long n, long prec)
{
    long len1, len2;

    len1 = poly1->length;
    len2 = poly2->length;

    if (len1 == 0 || len2 == 0 || n == 0)
    {
        fmpcb_poly_zero(res);
        return;
    }

    n = FLINT_MIN((len1 + len2 - 1), n);
    len1 = FLINT_MIN(len1, n);
    len2 = FLINT_MIN(len2, n);

    if (res == poly1 || res == poly2)
    {
        fmpcb_poly_t t;
        fmpcb_poly_init2(t, n);
        _fmpcb_poly_mullow_transpose(t->coeffs, poly1->coeffs, len1,
                                poly2->coeffs, len2, n, prec);
        fmpcb_poly_swap(res, t);
        fmpcb_poly_clear(t);
    }
    else
    {
        fmpcb_poly_fit_length(res, n);
        _fmpcb_poly_mullow_transpose(res->coeffs, poly1->coeffs, len1,
                                poly2->coeffs, len2, n, prec);
    }

    _fmpcb_poly_set_length(res, n);
    _fmpcb_poly_normalise(res);
}

