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

#include "arb.h"

static void
_arf_atan(arf_t z, const arf_t x, long prec, arf_rnd_t rnd)
{
    mpfr_t xf, zf;
    mp_ptr zptr, tmp;
    mp_srcptr xptr;
    mp_size_t xn, zn, val;
    TMP_INIT;
    TMP_START;

    zn = (prec + FLINT_BITS - 1) / FLINT_BITS;
    tmp = TMP_ALLOC(zn * sizeof(mp_limb_t));

    ARF_GET_MPN_READONLY(xptr, xn, x);

    xf->_mpfr_d = (mp_ptr) xptr;
    xf->_mpfr_prec = xn * FLINT_BITS;
    xf->_mpfr_sign = ARF_SGNBIT(x) ? -1 : 1;
    xf->_mpfr_exp = ARF_EXP(x);

    zf->_mpfr_d = tmp;
    zf->_mpfr_prec = prec;
    zf->_mpfr_sign = 1;
    zf->_mpfr_exp = 0;

    mpfr_set_emin(MPFR_EMIN_MIN);
    mpfr_set_emax(MPFR_EMAX_MAX);

    mpfr_atan(zf, xf, arf_rnd_to_mpfr(rnd));

    val = 0;
    while (tmp[val] == 0)
        val++;

    ARF_GET_MPN_WRITE(zptr, zn - val, z);
    flint_mpn_copyi(zptr, tmp + val, zn - val);
    if (zf->_mpfr_sign < 0)
        ARF_NEG(z);

    fmpz_set_si(ARF_EXPREF(z), zf->_mpfr_exp);

    TMP_END;
}

void
arb_atan_arf(arb_t z, const arf_t x, long prec)
{
    if (arf_is_special(x))
    {
        if (arf_is_zero(x))
        {
            arb_zero(z);
        }
        else if (arf_is_nan(x))
        {
            arb_indeterminate(z);
        }
        else if (arf_is_pos_inf(x))
        {
            arb_const_pi(z, prec);
            arb_mul_2exp_si(z, z, -1);
        }
        else if (arf_is_neg_inf(x))
        {
            arb_const_pi(z, prec);
            arb_mul_2exp_si(z, z, -1);
            arb_neg(z, z);
        }
    }
    else
    {
        /* 2^(mag-1) <= |x| < 2^mag, mag = ARF_EXP(x) */
        fmpz_t mag;

        if (ARF_EXP(x) >= -(prec/3) - 2 && ARF_EXP(x) <= prec + 2)
        {
            _arf_atan(arb_midref(z), x, prec, ARB_RND);
            arf_mag_set_ulp(arb_radref(z), arb_midref(z), prec);
        }
        else if (fmpz_sgn(ARF_EXPREF(x)) < 0)
        {
            /* atan(x) = x + eps, |eps| < x^3 */
            fmpz_init(mag);
            fmpz_mul_ui(mag, ARF_EXPREF(x), 3);
            arb_set_arf(z, x);
            arb_set_round(z, z, prec);
            arb_add_error_2exp_fmpz(z, mag);
            fmpz_clear(mag);
        }
        else
        {
            /* atan(x) = pi/2 - eps, eps < 1/x <= 2^(1-mag) */
            /* TODO: also use atan(x) = pi/2 - 1/x + eps, eps < 1/x^3 */
            fmpz_init(mag);

            fmpz_neg(mag, ARF_EXPREF(x));
            fmpz_add_ui(mag, mag, 1);

            if (arf_sgn(x) > 0)
            {
                arb_const_pi(z, prec);
            }
            else
            {
                arb_const_pi(z, prec);
                arb_neg(z, z);
            }

            arb_mul_2exp_si(z, z, -1);
            arb_add_error_2exp_fmpz(z, mag);

            fmpz_clear(mag);
        }
    }
}

void
arb_atan(arb_t z, const arb_t x, long prec)
{
    if (arb_is_exact(x))
    {
        arb_atan_arf(z, arb_midref(x), prec);
    }
    else
    {
        mag_t t, u;

        mag_init(t);
        mag_init(u);

        arb_get_mag_lower(t, x);

        if (mag_is_zero(t))
        {
            mag_set(t, arb_radref(x));
        }
        else
        {
            mag_mul_lower(t, t, t);
            mag_one(u);
            mag_add_lower(t, t, u);
            mag_div(t, arb_radref(x), t);
        }

        arb_atan_arf(z, arb_midref(x), prec);
        mag_add(arb_radref(z), arb_radref(z), t);

        mag_clear(t);
        mag_clear(u);
    }
}

