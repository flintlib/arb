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
_fmprb_const_log_sqrt2pi(fmprb_t t, long prec)
{
    fmprb_const_pi(t, prec + 2);
    fmprb_mul_2exp_si(t, t, 1);
    fmprb_log(t, t, prec);
    fmprb_mul_2exp_si(t, t, -1);
}

DEF_CACHED_CONSTANT(fmprb_const_log_sqrt2pi, _fmprb_const_log_sqrt2pi)
