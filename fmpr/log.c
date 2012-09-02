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

static mpfr_rnd_t fmpr_rnd_to_mpfr(fmpr_rnd_t rnd)
{
    if (rnd == FMPR_RND_NEAR) return MPFR_RNDN;
    if (rnd == FMPR_RND_FLOOR) return MPFR_RNDD;
    if (rnd == FMPR_RND_CEIL) return MPFR_RNDU;
    if (rnd == FMPR_RND_DOWN) return MPFR_RNDZ;
    if (rnd == FMPR_RND_UP) return MPFR_RNDA;
    abort();
}

long
fmpr_log(fmpr_t y, const fmpr_t x, long prec, fmpr_rnd_t rnd)
{
    if (fmpr_is_special(x))
    {
        if (fmpr_is_zero(x))
            fmpr_neg_inf(y);
        else if (fmpr_is_pos_inf(x))
            fmpr_pos_inf(y);
        else
            fmpr_nan(y);

        return FMPR_RESULT_EXACT;
    }

    if (fmpr_sgn(x) < 0)
    {
        fmpr_nan(y);
        return FMPR_RESULT_EXACT;
    }

    if (fmpr_is_one(x))
    {
        fmpr_zero(y);
        return FMPR_RESULT_EXACT;
    }

    {
        mpfr_t t, u;
        long r;

        mpfr_init2(t, 1 + fmpz_bits(fmpr_manref(x)));
        mpfr_init2(u, FLINT_MAX(2, prec));

        fmpr_get_mpfr(t, x, MPFR_RNDD);

        mpfr_log(u, t, fmpr_rnd_to_mpfr(rnd));

        fmpr_set_mpfr(y, u);

        r = prec - fmpz_bits(fmpr_manref(y));

        mpfr_clear(t);
        mpfr_clear(u);

        return r;
    }
}
