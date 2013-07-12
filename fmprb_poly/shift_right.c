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

    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "fmprb_poly.h"

void
_fmprb_poly_shift_right(fmprb_struct * res, const fmprb_struct * poly, long len, long n)
{
    long i;

    /* Copy in forward order to avoid writing over unshifted coefficients */
    if (res != poly)
    {
        for (i = 0; i < len - n; i++)
            fmprb_set(res + i, poly + n + i);
    }
    else
    {
        for (i = 0; i < len - n; i++)
            fmprb_swap(res + i, res + n + i);
    }

}

void
fmprb_poly_shift_right(fmprb_poly_t res, const fmprb_poly_t poly, long n)
{
    if (n == 0)
    {
        fmprb_poly_set(res, poly);
        return;
    }

    if (poly->length <= n)
    {
        fmprb_poly_zero(res);
        return;
    }

    fmprb_poly_fit_length(res, poly->length - n);
    _fmprb_poly_shift_right(res->coeffs, poly->coeffs, poly->length, n);
    _fmprb_poly_set_length(res, poly->length - n);
}

