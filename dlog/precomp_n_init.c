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
#include "math.h"

/* group of order n modulo mod, mod a prime and no information on n */
void
dlog_precomp_n_init(dlog_precomp_t pre, ulong a, ulong mod, ulong n, ulong num)
{
    if (n % 2 && n_is_probabprime(n))
        dlog_precomp_p_init(pre, a, mod, n, num);
    else {
        if (n < DLOG_TABLE_N_LIM)
        {
           dlog_precomp_small_init(pre, a, mod, n, num);
        }
        else
        {
            if (n < DLOG_BSGS_LIM)
            {
                ulong m;
                m = dlog_bsgs_size(n, num);
                pre->type = DLOG_BSGS;
                pre->cost = dlog_bsgs_init(pre->t.bsgs, a, mod, n, m);
            } else {
                pre->type = DLOG_CRT;
                pre->cost = dlog_crt_init(pre->t.crt, a, mod, n, num);
            }
        }
    }
}
