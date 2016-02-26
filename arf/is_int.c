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

    Copyright (C) 2016 Fredrik Johansson

******************************************************************************/

#include "arf.h"

int
arf_is_int(const arf_t x)
{
    mp_size_t xn;
    mp_srcptr xp;
    slong exp, c;

    exp = ARF_EXP(x);

    if (ARF_IS_SPECIAL(x))
        return exp == ARF_EXP_ZERO;

    if (COEFF_IS_MPZ(exp))
        return mpz_sgn(COEFF_TO_PTR(exp)) > 0;

    ARF_GET_MPN_READONLY(xp, xn, x);
    count_trailing_zeros(c, xp[0]);
    return exp - xn * FLINT_BITS + c >= 0;
}

