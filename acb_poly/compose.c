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

#include "acb_poly.h"

void
_acb_poly_compose(acb_ptr res,
    acb_srcptr poly1, long len1,
    acb_srcptr poly2, long len2, long prec)
{
    if (len2 == 1)
    {
        _acb_poly_evaluate(res, poly1, len1, poly2, prec);
    }
    else if (_acb_vec_is_zero(poly2, len2 - 1)) /* poly2 is a monomial */
    {
        long i;
        acb_t t;

        acb_init(t);
        acb_set(t, poly2 + len2 - 1);
        acb_set_round(res, poly1, prec);

        for (i = 1; i < len1; i++)
        {
            acb_mul(res + i * (len2 - 1), poly1 + i, t, prec);
            if (i + 1 < len1)
                acb_mul(t, t, poly2 + len2 - 1, prec);

            _acb_vec_zero(res + (i - 1) * (len2 - 1) + 1, len2 - 2);
        }

        acb_clear(t);
    }
    else if (len1 <= 7)
    {
        _acb_poly_compose_horner(res, poly1, len1, poly2, len2, prec);
    }
    else
    {
        _acb_poly_compose_divconquer(res, poly1, len1, poly2, len2, prec);
    }
}

void acb_poly_compose(acb_poly_t res,
              const acb_poly_t poly1, const acb_poly_t poly2, long prec)
{
    const long len1 = poly1->length;
    const long len2 = poly2->length;
    
    if (len1 == 0)
    {
        acb_poly_zero(res);
    }
    else if (len1 == 1 || len2 == 0)
    {
        acb_poly_set_acb(res, poly1->coeffs);
    }
    else
    {
        const long lenr = (len1 - 1) * (len2 - 1) + 1;
        
        if (res != poly1 && res != poly2)
        {
            acb_poly_fit_length(res, lenr);
            _acb_poly_compose(res->coeffs, poly1->coeffs, len1,
                                                   poly2->coeffs, len2, prec);
        }
        else
        {
            acb_poly_t t;
            acb_poly_init2(t, lenr);
            _acb_poly_compose(t->coeffs, poly1->coeffs, len1,
                                                 poly2->coeffs, len2, prec);
            acb_poly_swap(res, t);
            acb_poly_clear(t);
        }

        _acb_poly_set_length(res, lenr);
        _acb_poly_normalise(res);
    }
}
