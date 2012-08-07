/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "ufloat.h"

/* upper bound for log(2) */
#define LOG2 2977044472UL
#define LOG2_PREC 32

void
ufloat_log(ufloat_t z, const ufloat_t x)
{
    if (ufloat_cmp_one(x) <= 0)
        ufloat_zero(z);
    else
    {
        mp_limb_t hi, lo;

        umul_ppmm(hi, lo, x->exp, LOG2);
        ufloat_set_ll_2exp(z, hi, lo, -LOG2_PREC);
    }
}
