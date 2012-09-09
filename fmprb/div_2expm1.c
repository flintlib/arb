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

/* sets y = x / (2^n - 1) */
/* XXX: for huge n, use better algorithm / fix overflow in the loop */
void
fmprb_div_2expm1_ui(fmprb_t y, const fmprb_t x, ulong n, long prec)
{
    fmpz_t t;
    fmprb_t u;
    long i;

    fmpz_init(t);
    fmprb_init(u);

    for (i = prec; i >= 0; i -= n)
        fmpz_setbit(t, i);

    fmpr_set_fmpz(fmprb_midref(u), t);
    fmpr_one(fmprb_radref(u));

    fmprb_mul(y, x, u, prec);
    fmprb_mul_2exp_si(y, y, -prec-n);

    fmpz_clear(t);
    fmprb_clear(u);
}
