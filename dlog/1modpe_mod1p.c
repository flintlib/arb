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

    Copyright (C) 2016 Pascal Molin

******************************************************************************/

#include "dlog.h"

ulong
dlog_1modpe_mod1p(ulong b1, ulong p, ulong e, ulong inv1p, nmod_t pe)
{
    int f;
    ulong x, xf, pf, pf1;
    pf1 = 1;
    pf = p;
    x = 0;
    for (f = 1; f < e; f++)
    {      
        if (b1 % pf != 1)
            abort();
        xf = (b1 - 1) / pf;
        xf = (xf % p) * pf1;
        x += xf;
        b1 = nmod_mul(b1, nmod_pow_ui(inv1p, xf, pe), pe);
        pf1 = pf;
        pf *= p;
    }
    return x;
}
