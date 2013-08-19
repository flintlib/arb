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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include <math.h>
#include "fmprb_poly.h"

void
_fmprb_poly_binomial_transform_convolution(fmprb_ptr b, fmprb_srcptr a, long alen, long len, long prec)
{
    long i, s;
    fmprb_ptr c, d;

    alen = FLINT_MIN(alen, len);

    c = _fmprb_vec_init(alen);
    d = _fmprb_vec_init(len);

    /* Hack: rescale x -> 2^s x such that (2^s)^n / n! ~= 1,
       improving efficiency of block multiplication. The proper
       solution is to add automatic near-optimal scaling in the
       block multiplication. */
    s = (log(len) - 1) * 1.44269504088896 + 0.5;

    _fmprb_poly_borel_transform(c, a, alen, prec);
    for (i = 1; i < alen; i += 2)
        fmprb_neg(c + i, c + i);

    fmprb_one(d);
    for (i = 1; i < len; i++)
        fmprb_div_ui(d + i, d + i - 1, i, prec);

    for (i = 0; i < alen; i++)
        fmprb_mul_2exp_si(c + i, c + i, i * s);
    for (i = 0; i < len; i++)
        fmprb_mul_2exp_si(d + i, d + i, i * s);

    _fmprb_poly_mullow(b, d, len, c, alen, len, prec);

    for (i = 0; i < len; i++)
        fmprb_mul_2exp_si(b + i, b + i, -i * s);

    _fmprb_poly_inv_borel_transform(b, b, len, prec);

    _fmprb_vec_clear(c, alen);
    _fmprb_vec_clear(d, len);
}

void
fmprb_poly_binomial_transform_convolution(fmprb_poly_t b, const fmprb_poly_t a, long len, long prec)
{
    if (len == 0 || a->length == 0)
    {
        fmprb_poly_zero(b);
        return;
    }

    if (b == a)
    {
        fmprb_poly_t c;
        fmprb_poly_init2(c, len);
        _fmprb_poly_binomial_transform_convolution(c->coeffs, a->coeffs, a->length, len, prec);
        fmprb_poly_swap(b, c);
        fmprb_poly_clear(c);
    }
    else
    {
        fmprb_poly_fit_length(b, len);
        _fmprb_poly_binomial_transform_convolution(b->coeffs, a->coeffs, a->length, len, prec);
    }

    _fmprb_poly_set_length(b, len);
    _fmprb_poly_normalise(b);
}

