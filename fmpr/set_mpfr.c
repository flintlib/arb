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

void
fmpr_set_mpfr(fmpr_t x, const mpfr_t y)
{
    if (!mpfr_regular_p(y))
    {
        if (mpfr_zero_p(y)) fmpr_zero(x);
        else if (mpfr_inf_p(y) && mpfr_sgn(y) > 0) fmpr_pos_inf(x);
        else if (mpfr_inf_p(y) && mpfr_sgn(y) < 0) fmpr_neg_inf(x);
        else fmpr_nan(x);
    }
    else
    {
        __mpz_struct * z = _fmpz_promote(fmpr_manref(x));
        fmpz_set_si(fmpr_expref(x), mpfr_get_z_2exp(z, y));
        _fmpz_demote_val(fmpr_manref(x));
        _fmpr_normalise(fmpr_manref(x), fmpr_expref(x), y->_mpfr_prec, FMPR_RND_DOWN);
    }
}

