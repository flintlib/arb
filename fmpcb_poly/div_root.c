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

#include "fmpcb_poly.h"

void
_fmpcb_poly_div_root(fmpcb_struct * Q, fmpcb_t R, const fmpcb_struct * A,
    long len, const fmpcb_t c, long prec)
{
    fmpcb_t r, t;
    long i;

    if (len < 2)
    {
        fmpcb_zero(R);
        return;
    }

    fmpcb_init(r);
    fmpcb_init(t);

    fmpcb_set(t, A + len - 2);
    fmpcb_set(Q + len - 2, A + len - 1);
    fmpcb_set(r, Q + len - 2);

    /* TODO: avoid the extra assignments (but still support aliasing)  */
    for (i = len - 2; i > 0; i--)
    {
        fmpcb_mul(r, r, c, prec);
        fmpcb_add(r, r, t, prec);
        fmpcb_set(t, A + i - 1);
        fmpcb_set(Q + i - 1, r);
    }

    fmpcb_mul(r, r, c, prec);
    fmpcb_add(R, r, t, prec);
}
