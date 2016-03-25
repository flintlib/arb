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

void
dlog_1modpe_rec_init(dlog_1modpe_rec_t t, ulong a1, ulong p, ulong e, nmod_t pe)
{
    if (e == 1)
    {
        t->inv1p = 1;
        t->invloga1 = 0;
    }
    else
    {
        ulong loga1;
        if (a1 == 1)
            abort();
        t->inv1p = nmod_inv(1 + p, pe); /* 1 - p + p^2 - ... */
        loga1 = dlog_1modpe_mod1p(a1, p, e, t->inv1p, pe);
        /* only need inverse mod p^(e-1) but does not hurt */
        t->invloga1 = nmod_inv(loga1, pe);
    }
}
