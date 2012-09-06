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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "fmprb.h"

/* TODO: add log2(e) bits to the working precision? */

void
fmprb_pow_ui(fmprb_t y, const fmprb_t b, ulong e, long prec)
{
    long i;

    if (y == b)
    {
        fmprb_t t;
        fmprb_init(t);
        fmprb_set(t, b);
        fmprb_pow_ui(y, t, e, prec);
        fmprb_clear(t);
        return;
    }

    if (e == 0)
    {
        fmprb_set_ui(y, 1UL);
        return;
    }

    fmprb_set(y, b);

    for (i = FLINT_BIT_COUNT(e) - 2; i >= 0; i--)
    {
        fmprb_mul(y, y, y, prec);
        if (e & (1UL<<i))
            fmprb_mul(y, y, b, prec);
    }
}

void
fmprb_ui_pow_ui(fmprb_t y, ulong b, ulong e, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_set_ui(t, b);
    fmprb_pow_ui(y, t, e, prec);
    fmprb_clear(t);
}

void
fmprb_si_pow_ui(fmprb_t y, long b, ulong e, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_set_si(t, b);
    fmprb_pow_ui(y, t, e, prec);
    fmprb_clear(t);
}
