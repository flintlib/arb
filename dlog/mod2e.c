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
dlog_mod2e(const dlog_modpe_t t, ulong b1)
{
    if (t->e == 2)
        return (b1 % 4) == 3;
    else
    {
        slong f;
        ulong pf1, pf, x, xf;
        pf1 = 1;
        pf = 4;
        x = 0;
        for (f = 2; f < t->e; f++)
        {
            if (b1 % pf != 1)
            {
                flint_printf("ERROR dlog_mod2e: %wu %% %wu != 1 mod %wu\n\n",
                        b1, pf, t->pe.n);
                abort();
            }
            xf = (b1 - 1) / pf;
            xf = (f == 2) ? xf % 4 : (xf % 2) * (pf1 / 2);
            b1 = nmod_mul(b1, nmod_pow_ui(t->inva, xf, t->pe), t->pe);
            x += xf;
            pf1 = pf;
            pf *= 2;
        }
        return x;
    }
}
