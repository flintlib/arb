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

void
fmprb_addmul(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_mul(t, x, y, prec);
    fmprb_add(z, z, t, prec);
    fmprb_clear(t);
}

void
fmprb_addmul_ui(fmprb_t z, const fmprb_t x, ulong y, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_mul_ui(t, x, y, prec);
    fmprb_add(z, z, t, prec);
    fmprb_clear(t);
}

void
fmprb_addmul_si(fmprb_t z, const fmprb_t x, long y, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_mul_si(t, x, y, prec);
    fmprb_add(z, z, t, prec);
    fmprb_clear(t);
}

void
fmprb_addmul_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_mul_fmpz(t, x, y, prec);
    fmprb_add(z, z, t, prec);
    fmprb_clear(t);
}
