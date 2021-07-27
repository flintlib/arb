/*
    Copyright (C) 2021 Matthias Gessinger

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_graeffe_transform(arb_ptr b, arb_srcptr a, slong len, slong prec)
{
    if (len == 0)
        return;

    slong q, i;

    q = (len - 1) / 2 + 1;
    arb_ptr pe = _arb_vec_init(q);
    arb_ptr po = _arb_vec_init(len);

    for (i = len - 1; i >= 0; i--)
    {
        if (i % 2 == 0)
            arb_set(pe + i / 2, a + i);
        else
            arb_set(po + i / 2, a + i);
    }

    _arb_poly_mul(b, po, q, po, q, prec);
    _arb_poly_shift_left(b, b, len - 1, 1);
    _arb_poly_mul(po, pe, q, pe, q, prec);
    _arb_vec_sub(b, po, b, len, prec);

    _arb_vec_clear(pe, q);
    _arb_vec_clear(po, len);

    if (len % 2 == 0)
	    _arb_vec_neg(b, b, len);
}

void
arb_poly_graeffe_transform(arb_poly_t b, const arb_poly_t a, slong prec)
{
    arb_poly_fit_length(b, a->length);
    _arb_poly_graeffe_transform(b->coeffs, a->coeffs, a->length, prec);
    _arb_poly_set_length(b, a->length);
}
