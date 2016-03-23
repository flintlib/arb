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
dlog_precomp_small_init(dlog_precomp_t pre, ulong a, ulong mod, ulong n, ulong num)
{
    if (n <= 3)
    {
        pre->type = DLOG_23;
        pre->cost = dlog_order23_init(pre->t.order23, a); 
    }
    else
    {
        if ( mod < DLOG_TABLE_LIM )
        {
            pre->type = DLOG_TABLE;
            pre->cost = dlog_table_init(pre->t.table, a, mod);
        }
        else
        {
            pre->type = DLOG_BSGS;
            pre->cost = dlog_bsgs_init(pre->t.bsgs, a, mod, n, n);
        }
    }
}
