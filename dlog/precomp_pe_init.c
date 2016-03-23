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
dlog_precomp_pe_init(dlog_precomp_t pre, ulong a, ulong mod, ulong p, ulong e, ulong pe, ulong num)
{
    if ( pe < DLOG_TABLE_PE_LIM )
    {
        dlog_precomp_small_init(pre, a, mod, pe, num);
    }
    else
    {
        if ( e == 1)
        {
            dlog_precomp_p_init(pre, a, mod, p, num);
        }
        else
        {
            pre->type = DLOG_POWER;
            pre->cost = dlog_power_init(pre->t.power, a, mod, p, e, num);
        }
    }
}
