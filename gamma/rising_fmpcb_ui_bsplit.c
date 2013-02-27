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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "gamma.h"

void
gamma_rising_fmpcb_ui_bsplit(fmpcb_t y, const fmpcb_t x, ulong n, long prec)
{
    if (n < 8)
    {
        gamma_rising_fmpcb_ui_bsplit_simple(y, x, n, prec);
    }
    else if (prec < 250 || n < 250000 / prec)
    {
        gamma_rising_fmpcb_ui_bsplit_eight(y, x, n, prec);
    }
    else
    {
        ulong step, s1, s2;

        /* experimental fit */
        s1 = 2 * sqrt(n);
        s2 = 10 * pow(prec, 0.25);
        step = FLINT_MIN(s1, s2);

        gamma_rising_fmpcb_ui_bsplit_rectangular(y, x, n, step, prec);
    }
}

