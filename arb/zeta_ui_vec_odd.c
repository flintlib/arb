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

    Copyright (C) 2012, 2013 Fredrik Johansson

******************************************************************************/

#include "arb.h"

void
arb_zeta_ui_vec_odd(arb_ptr x, ulong start, long num, long prec)
{
    long i, num_borwein;
    ulong cutoff;

    cutoff = 40 + 0.3 * prec;

    if (cutoff > start)
    {
        num_borwein = 1 + (cutoff - start) / 2;
        num_borwein = FLINT_MIN(num_borwein, num);
    }
    else
        num_borwein = 0;

    arb_zeta_ui_vec_borwein(x, start, num_borwein, 2, prec);
    for (i = num_borwein; i < num; i++)
        arb_zeta_ui(x + i, start + 2 * i, prec);
}

