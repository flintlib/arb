/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

#define CUTOFF 4

void
_acb_poly_mullow(acb_ptr res,
    acb_srcptr poly1, slong len1,
    acb_srcptr poly2, slong len2, slong n, slong prec)
{
    if (n < CUTOFF || len1 < CUTOFF || len2 < CUTOFF)
        _acb_poly_mullow_classical(res, poly1, len1, poly2, len2, n, prec);
    else
        _acb_poly_mullow_transpose(res, poly1, len1, poly2, len2, n, prec);
}

void
acb_poly_mullow(acb_poly_t res, const acb_poly_t poly1,
                                            const acb_poly_t poly2,
                                                slong n, slong prec)
{
    slong len1, len2;

    len1 = poly1->length;
    len2 = poly2->length;

    if (len1 == 0 || len2 == 0 || n == 0)
    {
        acb_poly_zero(res);
        return;
    }

    n = FLINT_MIN((len1 + len2 - 1), n);
    len1 = FLINT_MIN(len1, n);
    len2 = FLINT_MIN(len2, n);

    if (res == poly1 || res == poly2)
    {
        acb_poly_t t;
        acb_poly_init2(t, n);
        _acb_poly_mullow(t->coeffs, poly1->coeffs, len1,
                                poly2->coeffs, len2, n, prec);
        acb_poly_swap(res, t);
        acb_poly_clear(t);
    }
    else
    {
        acb_poly_fit_length(res, n);
        _acb_poly_mullow(res->coeffs, poly1->coeffs, len1,
                                poly2->coeffs, len2, n, prec);
    }

    _acb_poly_set_length(res, n);
    _acb_poly_normalise(res);
}

