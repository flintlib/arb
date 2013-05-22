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

#include "elefun.h"

int
elefun_exp_precomp(fmprb_t z, const fmprb_t x, long prec, int minus_one)
{
    fmpz_t zfixed, zfixed_err, exponent, xfixed, xfixed_err;
    fmpr_t terr;
    long mag, fixed_wp;

    /* require radius < 0.25 */
    if (fmpr_cmp_2exp_si(fmprb_radref(x), -2) >= 0)
        return 0;

    if (fmpr_cmpabs_2exp_si(fmprb_midref(x), 128) >= 0)
        return 0;

    if (fmpr_cmpabs_2exp_si(fmprb_midref(x), -prec) < 0)
        return 0;

    fixed_wp = prec + 2 * FLINT_BIT_COUNT(prec);

    mag = *fmpr_expref(fmprb_midref(x)) + fmpz_bits(fmpr_manref(fmprb_midref(x)));

    if (mag > 0)
        fixed_wp += mag;   /* for argument reduction */
    else if (minus_one)
        fixed_wp += -mag;  /* cancellation */

    if (fixed_wp > EXP_CACHE_PREC)
        return 0;

    fmpz_init(zfixed);
    fmpz_init(zfixed_err);
    fmpz_init(exponent);
    fmpz_init(xfixed);
    fmpz_init(xfixed_err);
    fmpr_init(terr);

    fmpr_get_fmpz_fixed_si(xfixed, fmprb_midref(x), -fixed_wp);

    /* convert radius to fixed-point number; need to add 1 for
       error rounding midpoint, 1 for rounding the error itself */
    fmpr_get_fmpz_fixed_si(xfixed_err, fmprb_radref(x), -fixed_wp);
    fmpz_add_ui(xfixed_err, xfixed_err, 2);

    /* compute exp(x) as a fixed-point number */
    elefun_exp_fixed_precomp(zfixed, zfixed_err, exponent, xfixed, xfixed_err, fixed_wp);
    fmpz_sub_ui(exponent, exponent, fixed_wp);

    /* convert back to floating-point interval */
    if (minus_one)
    {
        fmprb_set_fmpz_2exp(z, zfixed, exponent);
        fmprb_sub_ui(z, z, 1, prec);
    }
    else
    {
        fmprb_set_round_fmpz_2exp(z, zfixed, exponent, prec);
    }

    /* add error */
    fmpr_set_round_fmpz_2exp(terr, zfixed_err, exponent, FMPRB_RAD_PREC, FMPR_RND_UP);
    fmpr_add(fmprb_radref(z), fmprb_radref(z), terr, FMPRB_RAD_PREC, FMPR_RND_UP);

    fmpz_clear(zfixed);
    fmpz_clear(zfixed_err);
    fmpz_clear(exponent);
    fmpz_clear(xfixed);
    fmpz_clear(xfixed_err);
    fmpr_clear(terr);

    return 1;
}

