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

#include "acb_poly.h"

void
_acb_poly_reverse(acb_ptr res, acb_srcptr poly, long len, long n)
{
    if (res == poly)
    {
        long i;

        for (i = 0; i < n / 2; i++)
        {
            acb_struct t = res[i];
            res[i] = res[n - 1 - i];
            res[n - 1 - i] = t;
        }

        for (i = 0; i < n - len; i++)
            acb_zero(res + i);
    }
    else
    {
        long i;

        for (i = 0; i < n - len; i++)
            acb_zero(res + i);

        for (i = 0; i < len; i++)
            acb_set(res + (n - len) + i, poly + (len - 1) - i);
    }
}
