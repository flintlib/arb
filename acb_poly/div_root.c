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

#include "acb_poly.h"

void
_acb_poly_div_root(acb_ptr Q, acb_t R, acb_srcptr A,
    long len, const acb_t c, long prec)
{
    acb_t r, t;
    long i;

    if (len < 2)
    {
        acb_zero(R);
        return;
    }

    acb_init(r);
    acb_init(t);

    acb_set(t, A + len - 2);
    acb_set(Q + len - 2, A + len - 1);
    acb_set(r, Q + len - 2);

    /* TODO: avoid the extra assignments (but still support aliasing)  */
    for (i = len - 2; i > 0; i--)
    {
        acb_mul(r, r, c, prec);
        acb_add(r, r, t, prec);
        acb_set(t, A + i - 1);
        acb_set(Q + i - 1, r);
    }

    acb_mul(r, r, c, prec);
    acb_add(R, r, t, prec);
}
