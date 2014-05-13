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

#include "arb_poly.h"

void
_arb_poly_div_root(arb_ptr Q, arb_t R, arb_srcptr A,
    long len, const arb_t c, long prec)
{
    arb_t r, t;
    long i;

    if (len < 2)
    {
        arb_zero(R);
        return;
    }

    arb_init(r);
    arb_init(t);

    arb_set(t, A + len - 2);
    arb_set(Q + len - 2, A + len - 1);
    arb_set(r, Q + len - 2);

    /* TODO: avoid the extra assignments (but still support aliasing)  */
    for (i = len - 2; i > 0; i--)
    {
        arb_mul(r, r, c, prec);
        arb_add(r, r, t, prec);
        arb_set(t, A + i - 1);
        arb_set(Q + i - 1, r);
    }

    arb_mul(r, r, c, prec);
    arb_add(R, r, t, prec);
}
