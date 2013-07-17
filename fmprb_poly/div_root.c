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

    Copyright (C) 2011 William Hart
    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "fmprb_poly.h"

void
_fmprb_poly_div_root(fmprb_ptr Q, fmprb_t R, fmprb_srcptr A,
    long len, const fmprb_t c, long prec)
{
    fmprb_t r, t;
    long i;

    if (len < 2)
    {
        fmprb_zero(R);
        return;
    }

    fmprb_init(r);
    fmprb_init(t);

    fmprb_set(t, A + len - 2);
    fmprb_set(Q + len - 2, A + len - 1);
    fmprb_set(r, Q + len - 2);

    /* TODO: avoid the extra assignments (but still support aliasing)  */
    for (i = len - 2; i > 0; i--)
    {
        fmprb_mul(r, r, c, prec);
        fmprb_add(r, r, t, prec);
        fmprb_set(t, A + i - 1);
        fmprb_set(Q + i - 1, r);
    }

    fmprb_mul(r, r, c, prec);
    fmprb_add(R, r, t, prec);
}
