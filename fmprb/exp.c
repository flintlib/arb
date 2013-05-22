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

    Copyright (C) 2012, 2013 Fredrik Johansson

******************************************************************************/

#include "fmprb.h"
#include "elefun.h"

void
fmpr_exp_ubound(fmpr_t y, const fmpr_t x, long prec, long maglim)
{
    fmpz_t mag;

    if (fmpr_is_special(x))
    {
        if (fmpr_is_pos_inf(x))
            fmpr_pos_inf(y);
        else if (fmpr_is_zero(x))
            fmpr_one(y);
        else if (fmpr_is_neg_inf(x))
            fmpr_zero(y);
        else
            fmpr_pos_inf(y);
        return;
    }

    fmpz_init(mag);

    /* 2^(mag-1) <= |x| < 2^mag */
    fmpz_add_ui(mag, fmpr_expref(x), fmpz_bits(fmpr_manref(x)));

    /* normal range -- ok to just call mpfr */
    if (!COEFF_IS_MPZ(*mag) && (FLINT_ABS(*mag) <= 24))
    {
        fmpr_exp(y, x, prec, FMPR_RND_UP);
    }
    /* close to zero */
    else if (fmpz_sgn(mag) < 0)
    {
        /* x < 0 ==> exp(x) < 1 */
        if (fmpz_sgn(fmpr_manref(x)) < 0)
        {
            fmpr_one(y);
        }
        /* x > 0 ==> exp(x) < 1 + 2x */
        else
        {
            fmpr_mul_2exp_si(y, x, 1);
            fmpr_add_ui(y, y, 1, prec, FMPR_RND_UP);
        }
    }
    /* overflow */
    else if (COEFF_IS_MPZ(*mag) || *mag > maglim)
    {
        if (fmpz_sgn(fmpr_manref(x)) > 0)
        {
            fmpr_pos_inf(y);
        }
        else
        {
            /*
               mag > maglim
               mag - 1 >= maglim

               x <= -2^(mag-1)
               x <= -2^maglim

               exp(x) < 2^(-2^(mag-1)) <= 2^(-2^maglim)
            */
            fmpz_set_si(fmpr_expref(y), -1);
            fmpz_mul_2exp(fmpr_expref(y), fmpr_expref(y), maglim);
            fmpz_one(fmpr_manref(y));
        }
    }
    /* +x or -x is huge */
    else
    {
        fmprb_t ln2, t, u;
        fmpz_t q;
        long wp;

        fmprb_init(ln2);
        fmprb_init(t);
        fmprb_init(u);
        fmpz_init(q);

        wp = prec + *mag + 10;

        fmprb_const_log2(ln2, wp);
        fmprb_set_fmpr(t, x);
        fmprb_div(u, t, ln2, wp);
        fmpr_get_fmpz(q, fmprb_midref(u), FMPR_RND_DOWN);
        fmprb_submul_fmpz(t, ln2, q, wp);

        fmpr_add(y, fmprb_midref(t), fmprb_radref(t), prec, FMPR_RND_CEIL);
        fmpr_exp(y, y, prec, FMPR_RND_UP);
        fmpr_mul_2exp_fmpz(y, y, q);

        fmprb_clear(ln2);
        fmprb_clear(t);
        fmprb_clear(u);
        fmpz_clear(q);
    }

    fmpz_clear(mag);
}

void
fmprb_exp_fmpr(fmprb_t z, const fmpr_t x, long prec, long maglim, int m1)
{
    long r, small_cutoff;
    fmpz_t mag;

    if (fmpr_is_special(x))
    {
        if (m1)
            fmpr_expm1(fmprb_midref(z), x, prec, FMPR_RND_DOWN);
        else
            fmpr_exp(fmprb_midref(z), x, prec, FMPR_RND_DOWN);
        fmpr_zero(fmprb_radref(z));
        return;
    }

    fmpz_init(mag);

    /* 2^(mag-1) <= |x| < 2^mag */
    fmpz_add_ui(mag, fmpr_expref(x), fmpz_bits(fmpr_manref(x)));

    /* magnitude for approximation by 1 + x */
    if (m1)
        small_cutoff = -FLINT_MAX(24, 4 + prec);
    else
        small_cutoff = -FLINT_MAX(24, 4 + prec / 2);

    /* standard case */
    if (!COEFF_IS_MPZ(*mag) && (*mag > small_cutoff) && (*mag <= 24))
    {
        if (m1)
            r = fmpr_expm1(fmprb_midref(z), x, prec, FMPR_RND_DOWN);
        else
            r = fmpr_exp(fmprb_midref(z), x, prec, FMPR_RND_DOWN);
        fmpr_set_error_result(fmprb_radref(z), fmprb_midref(z), r);
    }
    /* close to zero */
    else if (fmpz_sgn(mag) < 0)
    {
        /* exp(x) = 1 + x + eps,  eps < x^2 < 2^(2mag) */
        r = fmpr_add_ui(fmprb_midref(z), x, m1 ? 0 : 1, prec, FMPR_RND_DOWN);
        fmpr_set_error_result(fmprb_radref(z), fmprb_midref(z), r);
        fmpz_mul_2exp(mag, mag, 1);
        fmprb_add_error_2exp_fmpz(z, mag);
    }
    /* overflow */
    else if (COEFF_IS_MPZ(*mag) || *mag > maglim)
    {
        if (fmpz_sgn(fmpr_manref(x)) > 0)
        {
            fmpr_pos_inf(fmprb_midref(z));
            fmpr_pos_inf(fmprb_radref(z));
        }
        else
        {
            /*
               mag > maglim
               mag - 1 >= maglim

               x <= -2^(mag-1)
               x <= -2^maglim

               0 < exp(x) < 2^(-2^(mag-1)) <= 2^(-2^maglim)
            */
            fmpz_set_si(fmpr_expref(fmprb_midref(z)), -1);
            fmpz_mul_2exp(fmpr_expref(fmprb_midref(z)), fmpr_expref(fmprb_midref(z)), maglim);
            fmpz_one(fmpr_manref(fmprb_midref(z)));
            fmpr_set(fmprb_radref(z), fmprb_midref(z));
            if (m1)
                fmprb_sub_ui(z, z, 1, prec);
        }
    }
    /* huge */
    else
    {
        fmprb_t ln2, t, u;
        fmpz_t q;
        long wp;

        fmprb_init(ln2);
        fmprb_init(t);
        fmprb_init(u);
        fmpz_init(q);

        wp = prec + *mag + 10;

        fmprb_const_log2(ln2, wp);
        fmprb_set_fmpr(t, x);
        fmprb_div(u, t, ln2, wp);
        fmpr_get_fmpz(q, fmprb_midref(u), FMPR_RND_DOWN);
        fmprb_submul_fmpz(t, ln2, q, wp);

        fmprb_exp(z, t, prec);
        fmprb_mul_2exp_fmpz(z, z, q);

        if (m1)
            fmprb_sub_ui(z, z, 1, prec);

        fmprb_clear(ln2);
        fmprb_clear(t);
        fmprb_clear(u);
        fmpz_clear(q);
    }

    fmpz_clear(mag);
}

void
_fmprb_exp(fmprb_t z, const fmprb_t x, long prec, int m1)
{
    long maglim = FLINT_MAX(128, 2 * prec);

    if (elefun_exp_precomp(z, x, prec, m1))
        return;

    if (fmprb_is_exact(x))
    {
        fmprb_exp_fmpr(z, fmprb_midref(x), prec, maglim, m1);
    }
    else
    {
        /* exp(a+b) - exp(a) = exp(a) * (exp(b)-1) <= b * exp(a+b) */
        fmpr_t t;
        fmpr_init(t);
        fmpr_add(t, fmprb_midref(x), fmprb_radref(x), prec, FMPR_RND_CEIL);
        fmpr_exp_ubound(t, t, FMPRB_RAD_PREC, maglim);
        fmpr_mul(t, t, fmprb_radref(x), FMPRB_RAD_PREC, FMPR_RND_UP);
        fmprb_exp_fmpr(z, fmprb_midref(x), prec, maglim, m1);
        fmpr_add(fmprb_radref(z), fmprb_radref(z), t, FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_clear(t);
    }
    fmprb_adjust(z);
}

void
fmprb_exp(fmprb_t z, const fmprb_t x, long prec)
{
    _fmprb_exp(z, x, prec, 0);
}

void
fmprb_expm1(fmprb_t z, const fmprb_t x, long prec)
{
    _fmprb_exp(z, x, prec, 1);
}

