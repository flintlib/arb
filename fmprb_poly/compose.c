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

    Copyright (C) 2010 William Hart
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "fmprb_poly.h"

void
_fmprb_poly_compose(fmprb_ptr res,
    fmprb_srcptr poly1, long len1,
    fmprb_srcptr poly2, long len2, long prec)
{
    if (len2 == 1)
    {
        _fmprb_poly_evaluate(res, poly1, len1, poly2, prec);
    }
    else if (_fmprb_vec_is_zero(poly2, len2 - 1)) /* poly2 is a monomial */
    {
        long i;
        fmprb_t t;

        fmprb_init(t);
        fmprb_set(t, poly2 + len2 - 1);
        fmprb_set_round(res, poly1, prec);

        for (i = 1; i < len1; i++)
        {
            fmprb_mul(res + i * (len2 - 1), poly1 + i, t, prec);
            if (i + 1 < len1)
                fmprb_mul(t, t, poly2 + len2 - 1, prec);

            _fmprb_vec_zero(res + (i - 1) * (len2 - 1) + 1, len2 - 2);
        }

        fmprb_clear(t);
    }
    else if (len1 <= 7)
    {
        _fmprb_poly_compose_horner(res, poly1, len1, poly2, len2, prec);
    }
    else
    {
        _fmprb_poly_compose_divconquer(res, poly1, len1, poly2, len2, prec);
    }
}

void fmprb_poly_compose(fmprb_poly_t res,
              const fmprb_poly_t poly1, const fmprb_poly_t poly2, long prec)
{
    const long len1 = poly1->length;
    const long len2 = poly2->length;
    
    if (len1 == 0)
    {
        fmprb_poly_zero(res);
    }
    else if (len1 == 1 || len2 == 0)
    {
        fmprb_poly_set_fmprb(res, poly1->coeffs);
    }
    else
    {
        const long lenr = (len1 - 1) * (len2 - 1) + 1;
        
        if (res != poly1 && res != poly2)
        {
            fmprb_poly_fit_length(res, lenr);
            _fmprb_poly_compose(res->coeffs, poly1->coeffs, len1,
                                                   poly2->coeffs, len2, prec);
        }
        else
        {
            fmprb_poly_t t;
            fmprb_poly_init2(t, lenr);
            _fmprb_poly_compose(t->coeffs, poly1->coeffs, len1,
                                                 poly2->coeffs, len2, prec);
            fmprb_poly_swap(res, t);
            fmprb_poly_clear(t);
        }

        _fmprb_poly_set_length(res, lenr);
        _fmprb_poly_normalise(res);
    }
}
