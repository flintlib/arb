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

#include "gamma.h"

void
gamma_taylor_bound_ratio(fmpr_t r, long n)
{
    fmprb_t t;
    fmpq_t q;

    fmprb_init(t);
    fmpq_init(q);
    fmpq_set_si(q, -5, 6);

    fmprb_set_ui(t, n);
    fmprb_pow_fmpq(t, t, q, FMPRB_RAD_PREC);
    fmprb_mul_ui(t, t, 3, FMPRB_RAD_PREC);

    fmpr_add(r, fmprb_midref(t), fmprb_radref(t), FMPRB_RAD_PREC, FMPR_RND_UP);

    fmprb_clear(t);
    fmpq_clear(q);
}

