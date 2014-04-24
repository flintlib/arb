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

#include "arf.h"

void
arf_set_mpfr(arf_t x, const mpfr_t y)
{
    if (!mpfr_regular_p(y))
    {
        if (mpfr_zero_p(y))
            arf_zero(x);
        else if (mpfr_inf_p(y) && mpfr_sgn(y) > 0)
            arf_pos_inf(x);
        else if (mpfr_inf_p(y) && mpfr_sgn(y) < 0)
            arf_neg_inf(x);
        else
            arf_nan(x);
    }
    else
    {
        mp_size_t n = (y->_mpfr_prec + FLINT_BITS - 1) / FLINT_BITS;
        arf_set_mpn(x, y->_mpfr_d, n, y->_mpfr_sign < 0);
        fmpz_set_si(ARF_EXPREF(x), y->_mpfr_exp);
    }
}

