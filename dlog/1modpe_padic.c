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

/* assume b = 1 mod p, not checked */
ulong
dlog_1modpe_padic(const dlog_1modpe_padic_t t, ulong b1)
{
    padic_t px;
    fmpz_t ix;
    ulong ux;

    if (b1 == 1)
        return 0;

    padic_init(px);
    fmpz_init(ix);

    padic_set_ui(px, b1, t->ctx);
    flint_printf("set %wu -> ", b1);
    padic_print(px, t->ctx);

    padic_log(px, px, t->ctx);
    flint_printf("\n\n compute log -> ");
    padic_print(px, t->ctx);

    flint_printf("\n\n 1/log(a^(p-1)) -> ");
    padic_print(t->invlog, t->ctx);

    padic_mul(px, px, t->invlog, t->ctx);
    flint_printf("\n\n divide by log(a^(p-1)) -> ");
    padic_print(px, t->ctx);

    padic_get_fmpz(ix, px, t->ctx);
    ux =  fmpz_get_ui(ix);
    flint_printf("\n\nlog_p(%wu)/log_p(a) = %wu\n", b1, ux);
    padic_clear(px);
    fmpz_clear(ix);
    return ux;
}
