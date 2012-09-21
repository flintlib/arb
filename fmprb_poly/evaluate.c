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

#include "fmprb_poly.h"

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

    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "fmprb_poly.h"

void
_fmprb_poly_evaluate(fmprb_t res, const fmprb_struct * f, long len,
                           const fmprb_t a, long prec)
{
    if (len == 0)
    {
        fmprb_zero(res);
    }
    else if (len == 1 || fmprb_is_zero(a))
    {
        fmprb_set(res, f);
    }
    else
    {
        long i = len - 1;
        fmprb_t t;

        fmprb_init(t);
        fmprb_set(res, f + i);
        for (i = len - 2; i >= 0; i--)
        {
            fmprb_mul(t, res, a, prec);
            fmprb_add(res, f + i, t, prec);
        }
        fmprb_clear(t);
    }
}

void
fmprb_poly_evaluate(fmprb_t res, const fmprb_poly_t f, const fmprb_t a, long prec)
{
    if (res == a)
    {
        fmprb_t t;
        fmprb_init(t);
        _fmprb_poly_evaluate(t, f->coeffs, f->length, a, prec);
        fmprb_swap(res, t);
        fmprb_clear(t);
    }
    else
        _fmprb_poly_evaluate(res, f->coeffs, f->length, a, prec);
}
