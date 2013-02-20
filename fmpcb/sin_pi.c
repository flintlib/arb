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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "fmpcb.h"

void
fmpcb_sin_pi(fmpcb_t r, const fmpcb_t z, long prec)
{
#define a fmpcb_realref(z)
#define b fmpcb_imagref(z)

    fmprb_t sa, ca, sb, cb;

    fmprb_init(sa);
    fmprb_init(ca);
    fmprb_init(sb);
    fmprb_init(cb);

    fmprb_sin_cos_pi(sa, ca, a, prec);
    fmprb_const_pi(cb, prec);
    fmprb_mul(cb, cb, b, prec);
    fmprb_sinh_cosh(sb, cb, cb, prec);

    fmprb_mul(fmpcb_realref(r), sa, cb, prec);
    fmprb_mul(fmpcb_imagref(r), sb, ca, prec);

    fmprb_clear(sa);
    fmprb_clear(ca);
    fmprb_clear(sb);
    fmprb_clear(cb);

#undef a
#undef b
}

