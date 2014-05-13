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

#include "arb_poly.h"

void
_arb_poly_compose_horner(arb_ptr res,
    arb_srcptr poly1, long len1,
    arb_srcptr poly2, long len2, long prec)
{
    if (len1 == 1)
    {
        arb_set(res, poly1);
    }
    else if (len2 == 1)
    {
        _arb_poly_evaluate(res, poly1, len1, poly2, prec);
    }
    else if (len1 == 2)
    {
        _arb_vec_scalar_mul(res, poly2, len2, poly1 + 1, prec);
        arb_add(res, res, poly1, prec);
    }
    else
    {
        const long alloc = (len1 - 1) * (len2 - 1) + 1;
        long i = len1 - 1, lenr = len2;
        arb_ptr t, t1, t2;
        t = _arb_vec_init(alloc);

        if (len1 % 2 == 0)
        {
            t1 = res;
            t2 = t;
        }
        else
        {
            t1 = t;
            t2 = res;
        }

        /* Perform the first two steps as one,
            "res = a(m) * poly2 + a(m-1)". */
        {
            _arb_vec_scalar_mul(t1, poly2, len2, poly1 + i, prec);
            i--;
            arb_add(t1 + 0, t1 + 0, poly1 + i, prec);
        }
        while (i--)
        {
            _arb_poly_mul(t2, t1, lenr, poly2, len2, prec);
            lenr += len2 - 1;
            {
                void *t_ = t1;
                t1 = t2;
                t2 = t_;
            }
            arb_add(t1 + 0, t1 + 0, poly1 + i, prec);
        }
        _arb_vec_clear(t, alloc);
    }
}

void arb_poly_compose_horner(arb_poly_t res,
              const arb_poly_t poly1, const arb_poly_t poly2, long prec)
{
    const long len1 = poly1->length;
    const long len2 = poly2->length;
    
    if (len1 == 0)
    {
        arb_poly_zero(res);
    }
    else if (len1 == 1 || len2 == 0)
    {
        arb_poly_set_arb(res, poly1->coeffs);
    }
    else
    {
        const long lenr = (len1 - 1) * (len2 - 1) + 1;
        
        if (res != poly1 && res != poly2)
        {
            arb_poly_fit_length(res, lenr);
            _arb_poly_compose_horner(res->coeffs, poly1->coeffs, len1,
                                                   poly2->coeffs, len2, prec);
        }
        else
        {
            arb_poly_t t;
            arb_poly_init2(t, lenr);
            _arb_poly_compose_horner(t->coeffs, poly1->coeffs, len1,
                                                 poly2->coeffs, len2, prec);
            arb_poly_swap(res, t);
            arb_poly_clear(t);
        }

        _arb_poly_set_length(res, lenr);
        _arb_poly_normalise(res);
    }
}
