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

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#ifndef ARB_H
#define ARB_H

#include "fmprb.h"
#include "mag.h"
#include "arf.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    arf_struct mid;
    mag_struct rad;
}
arb_struct;

typedef arb_struct arb_t[1];
typedef arb_struct * arb_ptr;
typedef const arb_struct * arb_srcptr;

#define ARB_MIDREF(x) (&(x)->mid)
#define ARB_RADREF(x) (&(x)->rad)

#define ARB_IS_LAGOM(x) (ARF_IS_LAGOM(ARB_MIDREF(x)) && MAG_IS_LAGOM(ARB_RADREF(x)))


static __inline__ void
arb_init(arb_t x)
{
    arf_init(ARB_MIDREF(x));
    mag_init(ARB_RADREF(x));
}

static __inline__ void
arb_clear(arb_t x)
{
    arf_clear(ARB_MIDREF(x));
    mag_clear(ARB_RADREF(x));
}

static __inline__ void
arb_zero(arb_t x)
{
    arf_zero(ARB_MIDREF(x));
    mag_zero(ARB_RADREF(x));
}

static __inline__ void
arb_set_fmprb(arb_t x, const fmprb_t y)
{
    arf_set_fmpr(ARB_MIDREF(x), fmprb_midref(y));
    mag_set_fmpr(ARB_RADREF(x), fmprb_radref(y));
}

static __inline__ void
arb_get_fmprb(fmprb_t x, const arb_t y)
{
    arf_get_fmpr(fmprb_midref(x), ARB_MIDREF(y));
    mag_get_fmpr(fmprb_radref(x), ARB_RADREF(y));
}

static __inline__ void
arb_printd(const arb_t y, long d)
{
    fmprb_t t;
    fmprb_init(t);
    arb_get_fmprb(t, y);
    fmprb_printd(t, d);
    fmprb_clear(t);
}

void arb_mul(arb_t z, const arb_t x, const arb_t y, long prec);

void arb_addmul(arb_t z, const arb_t x, const arb_t y, long prec);

#ifdef __cplusplus
}
#endif

#endif

