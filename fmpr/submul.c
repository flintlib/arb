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

#include "fmpr.h"

long
fmpr_submul(fmpr_t z, const fmpr_t x, const fmpr_t y, long prec, fmpr_rnd_t rnd)
{
    fmpr_t t;
    long r;
    fmpr_init(t);
    fmpr_mul(t, x, y, LONG_MAX, FMPR_RND_DOWN);
    r = fmpr_sub(z, z, t, prec, rnd);
    fmpr_clear(t);
    return r;
}
