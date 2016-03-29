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
dlog_precomp_clear(dlog_precomp_t pre)
{
    if (pre == NULL)
        return;
    switch (pre->type)
    {
        case DLOG_MODPE:
            dlog_modpe_clear(pre->t.modpe);
            break;
        case DLOG_CRT:
            dlog_crt_clear(pre->t.crt);
            break;
        case DLOG_POWER:
            dlog_power_clear(pre->t.power);
            break;
        case DLOG_TABLE:
            dlog_table_clear(pre->t.table);
            break;
        case DLOG_BSGS:
            dlog_bsgs_clear(pre->t.bsgs);
            break;
        case DLOG_23:
            dlog_order23_clear(pre->t.order23);
            break;
        default:
            abort();
            break;
    }
}
