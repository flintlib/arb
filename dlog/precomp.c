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
dlog_precomp(const dlog_precomp_t pre, ulong b)
{
    if (b == 1)
        return 0;
    switch (pre->type)
    {
      case DLOG_MODPE:
        return dlog_modpe(pre->t.modpe, b);
        break;
      case DLOG_CRT:
        return dlog_crt(pre->t.crt, b);
        break;
      case DLOG_POWER:
        return dlog_power(pre->t.power, b);
        break;
      case DLOG_TABLE:
        return dlog_table(pre->t.table, b);
        break;
      case DLOG_BSGS:
        return dlog_bsgs(pre->t.bsgs, b);
        break;
      default:
        abort();
        break;
    }
}
