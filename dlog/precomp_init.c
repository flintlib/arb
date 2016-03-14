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

#include "math.h"
#include "dlog.h"

#define TABLE_P_LIM 50
#define TABLE_PE_LIM 50
#define TABLE_N_LIM 50
#define BSGS_LIM  500 


/* log mod p^e */
void
dlog_precomp_modpe_init(dlog_precomp_t pre, ulong a, ulong p, ulong e, ulong pe, ulong num)
{
    if ( pe < TABLE_PE_LIM )
    {
        pre->type = DLOG_TABLE;
        dlog_table_init(pre->t.table, a, pe);
        return;
    }
    if (e > 1)
    {
        pre->type = DLOG_MODPE;
        dlog_modpe_init(pre->t.modpe, a, p, e, pe, num);
    }
    else
    {
        dlog_precomp_n_init(pre, a, p, p - 1, num);
    }
}

/* group of order n modulo mod, mod a prime and no information on n */
void
dlog_precomp_n_init(dlog_precomp_t pre, ulong a, ulong mod, ulong n, ulong num)
{
    if (n%2 && n_is_probabprime(n))
        dlog_precomp_p_init(pre, a, mod, n, num);
    else {
        if ( mod < TABLE_N_LIM )
        {
           pre->type = DLOG_TABLE;
           dlog_table_init(pre->t.table, a, mod);
        } else {
            if (n < BSGS_LIM)
            {
                ulong m;
                m = (2 * num < n) ? ceil(sqrt((double) n * num)) : n;
                pre->type = DLOG_BSGS;
                dlog_bsgs_init(pre->t.bsgs, a, mod, n, m);
            } else {
                pre->type = DLOG_CRT;
                dlog_crt_init(pre->t.crt, a, mod, n, num);
            }
        }
    }
}

/* we known the order is prime */
void
dlog_precomp_p_init(dlog_precomp_t pre, ulong a, ulong mod, ulong p, ulong num)
{
    if ( p < TABLE_P_LIM )
    {
        pre->type = DLOG_TABLE;
        dlog_table_init(pre->t.table, a, mod);
    }
    else
    {
        ulong m;
        m = (2 * num < p) ? ceil(sqrt((double) p * num)) : p;
        pre->type = DLOG_BSGS;
        dlog_bsgs_init(pre->t.bsgs, a, mod, p, m);
    }
}

void
dlog_precomp_pe_init(dlog_precomp_t pre, ulong a, ulong mod, ulong p, ulong e, ulong pe, ulong num)
{
    if ( pe < TABLE_PE_LIM )
    {
        pre->type = DLOG_TABLE;
        dlog_table_init(pre->t.table, a, mod);
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
            dlog_power_init(pre->t.power, a, mod, p, e, num);
        }
    }
}

void
dlog_precomp_clear(dlog_precomp_t pre)
{
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
      default:
        abort();
        break;
    }
}
