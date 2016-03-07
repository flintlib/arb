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

#include "acb_dirichlet.h"
#include "padic.h"

typedef struct
{
    ulong p;
    padic_ctx_t ctx;     /* padic context */
    padic_t invlog;
}
dlog_1modpe_struct;

typedef dlog_1modpe_struct dlog_1modpe_t[1];

void
dlog_1modpe_init(dlog_1modpe_t t, ulong a1, ulong p, ulong e)
{
    fmpz_t tmp;
    t->p = p;
    fmpz_init(tmp);
    padic_init(t->invlog);

    fmpz_set_ui(tmp, p);
    padic_ctx_init(t->ctx , tmp , 0 , e, PADIC_TERSE);

    padic_set_ui(t->invlog, a1, t->ctx);
    padic_inv(t->invlog, t->invlog, t->ctx);

    fmpz_clear(tmp);
}

void
dlog_1modpe_clear(dlog_1modpe_t t)
{
    padic_clear(t->invlog);
    padic_ctx_clear(t->ctx);
}

/* assume b = 1 mod p, not checked */
ulong
dlog_1modpe(const dlog_1modpe_t t, ulong b)
{
    padic_t px;
    fmpz_t ix;
    ulong ux;
    padic_init(px);
    fmpz_init(ix);

    padic_set_ui(px, b, t->ctx);
    padic_log(px, px, t->ctx);
    padic_mul(px, px, t->invlog, t->ctx);
    padic_get_fmpz(ix, px, t->ctx);
    ux =  fmpz_get_ui(ix);

    padic_clear(px);
    fmpz_clear(ix);
    return ux;
}
