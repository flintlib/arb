/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_mullow(arb_ptr res,
    arb_srcptr poly1, slong len1,
    arb_srcptr poly2, slong len2, slong n, slong prec)
{
    if (n == 1)
    {
        arb_mul(res, poly1, poly2, prec);
    }
    else if (n <= 8 || len1 <= 8 || len2 <= 8)
    {
        _arb_poly_mullow_classical(res, poly1, len1, poly2, len2, n, prec);
    }
    else
    {
        slong m, cutoff;

        len1 = FLINT_MIN(len1, n);
        len2 = FLINT_MIN(len2, n);
        m = FLINT_MAX(len1, len2);
        m = FLINT_MAX(m, n);

        if      (prec <= 128)   cutoff = 100;
        else if (prec <= 192)   cutoff = 55;
        else if (prec <= 256)   cutoff = 70;
        else if (prec <= 512)   cutoff = 50;
        else if (prec <= 1024)  cutoff = 35;
        else if (prec <= 2048)  cutoff = 25;
        else if (prec <= 4096)  cutoff = 35;
        else if (prec <= 8192)  cutoff = 25;
        else if (prec <= 16384) cutoff = 20;
        else                    cutoff = 15;

        if (poly1 == poly2 && prec >= 256)
            cutoff *= 1.25;

        if (m < cutoff)
            _arb_poly_mullow_classical(res, poly1, len1, poly2, len2, n, prec);
        else
            _arb_poly_mullow_block(res, poly1, len1, poly2, len2, n, prec);
    }
}

void
arb_poly_mullow(arb_poly_t res, const arb_poly_t poly1,
                        const arb_poly_t poly2, slong n, slong prec)
{
    slong len1, len2;

    len1 = poly1->length;
    len2 = poly2->length;

    if (len1 == 0 || len2 == 0 || n == 0)
    {
        arb_poly_zero(res);
        return;
    }

    n = FLINT_MIN((len1 + len2 - 1), n);
    len1 = FLINT_MIN(len1, n);
    len2 = FLINT_MIN(len2, n);

    /* Hack to avoid temporary allocations with first derivatives. */
    if (n <= 2 && !(len1 == 2 && len2 == 2))
    {
        arb_poly_fit_length(res, n);

        if (n == 1)
        {
            arb_mul(res->coeffs, poly1->coeffs, poly2->coeffs, prec);
        }
        else if (len2 == 1)
        {
            arb_mul(res->coeffs + 1, poly1->coeffs + 1, poly2->coeffs, prec);
            arb_mul(res->coeffs, poly1->coeffs, poly2->coeffs, prec);
        }
        else if (len1 == 1)
        {
            arb_mul(res->coeffs + 1, poly2->coeffs + 1, poly1->coeffs, prec);
            arb_mul(res->coeffs, poly2->coeffs, poly1->coeffs, prec);
        }
        else
        {
            abort();

            if (res == poly1 || res == poly2)
            {
                arb_t t;
                arb_init(t);
                arb_mul(t, poly1->coeffs, poly2->coeffs + 1, prec);
                arb_addmul(t, poly2->coeffs, poly1->coeffs + 1, prec);
                arb_mul(res->coeffs, poly1->coeffs, poly2->coeffs, prec);
                arb_swap(t, res->coeffs + 1);
                arb_clear(t);
            }
            else
            {
                arb_mul(res->coeffs, poly1->coeffs, poly2->coeffs, prec);
                arb_mul(res->coeffs + 1, poly1->coeffs, poly2->coeffs + 1, prec);
                arb_addmul(res->coeffs + 1, poly2->coeffs, poly1->coeffs + 1, prec);
            }
        }

        _arb_poly_set_length(res, n);
        _arb_poly_normalise(res);
        return;
    }

    if (res == poly1 || res == poly2)
    {
        arb_poly_t t;
        arb_poly_init2(t, n);
        _arb_poly_mullow(t->coeffs, poly1->coeffs, len1,
                                poly2->coeffs, len2, n, prec);
        arb_poly_swap(res, t);
        arb_poly_clear(t);
    }
    else
    {
        arb_poly_fit_length(res, n);
        _arb_poly_mullow(res->coeffs, poly1->coeffs, len1,
                                poly2->coeffs, len2, n, prec);
    }

    _arb_poly_set_length(res, n);
    _arb_poly_normalise(res);
}

