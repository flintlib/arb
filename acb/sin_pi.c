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

#include "acb.h"

void
acb_sin_pi(acb_t r, const acb_t z, long prec)
{
#define a acb_realref(z)
#define b acb_imagref(z)

    arb_t sa, ca, sb, cb;

    arb_init(sa);
    arb_init(ca);
    arb_init(sb);
    arb_init(cb);

    arb_sin_cos_pi(sa, ca, a, prec);
    arb_const_pi(cb, prec);
    arb_mul(cb, cb, b, prec);
    arb_sinh_cosh(sb, cb, cb, prec);

    arb_mul(acb_realref(r), sa, cb, prec);
    arb_mul(acb_imagref(r), sb, ca, prec);

    arb_clear(sa);
    arb_clear(ca);
    arb_clear(sb);
    arb_clear(cb);

#undef a
#undef b
}

