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

/* TODO: make  threadsafe; round output */

long fmprb_const_pi_cached_prec = 0;
fmprb_t fmprb_const_pi_cache;

void
fmprb_const_pi(fmprb_t x, long prec)
{
    if (fmprb_const_pi_cached_prec < prec)
    {
        if (fmprb_const_pi_cached_prec == 0)
            fmprb_init(fmprb_const_pi_cache);

        fmprb_const_pi_chudnovsky(fmprb_const_pi_cache, prec);
        fmprb_const_pi_cached_prec = prec;
    }

    fmprb_set(x, fmprb_const_pi_cache);
}
