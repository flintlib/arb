/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
_acb_poly_graeffe_transform(acb_ptr b, acb_srcptr a, slong len, slong prec)
{
    if (len == 0)
        return;

    slong q, i;

    q = (len-1)/2+1;
    acb_ptr pe = _acb_vec_init(q);
    acb_ptr po = _acb_vec_init(len);

    for (i = len-1; i >= 0; i--)
    {
        if (i % 2 == 0)
            acb_set(pe+i/2, a+i);
        else
            acb_set(po+i/2, a+i);
    }

    _acb_poly_mul(b, po, q, po, q, prec);
    _acb_poly_shift_left(b, b, len-1, 1);
    _acb_poly_mul(po, pe, q, pe, q, prec);
    _acb_vec_sub(b, po, b, len, prec);

    _acb_vec_clear(pe, q);
    _acb_vec_clear(po, len);

    if (len % 2 == 0)
	    _acb_vec_neg(b, b, len);
}

void
acb_poly_graeffe_transform(acb_poly_t b, const acb_poly_t a, slong prec)
{
    acb_poly_fit_length(b, a->length);
    _acb_poly_graeffe_transform(b->coeffs, a->coeffs, a->length, prec);
    _acb_poly_set_length(b, a->length);
}
