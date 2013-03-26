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

void
_fmprb_pow_exp(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_log(t, x, prec);
    fmprb_mul(t, t, y, prec);
    fmprb_exp(z, t, prec);
    fmprb_clear(t);
}

void
fmprb_pow(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)
{
    if (fmprb_is_exact(y))
    {
        if (fmpr_is_zero(fmprb_midref(y)))
        {
            fmprb_one(z);
            return;
        }

        if (!fmpr_is_special(fmprb_midref(y)) && !fmpr_is_special(fmprb_midref(x)))
        {
            const fmpz * exp_exp = fmpr_expref(fmprb_midref(y));

            /* smallish integer powers and square roots */
            if (!COEFF_IS_MPZ(*exp_exp) && (*exp_exp >= -1L))
            {
                long exp_bits;

                exp_bits = *exp_exp + fmpz_bits(fmpr_manref(fmprb_midref(y)));

                if (exp_bits < 64)
                {
                    fmpz_t e;
                    fmpz_init(e);

                    if (*exp_exp == -1L)
                    {
                        fmprb_sqrt(z, x, prec + exp_bits);
                        fmpz_set(e, fmpr_manref(fmprb_midref(y)));
                        fmprb_pow_fmpz_binexp(z, z, e, prec);
                    }
                    else
                    {
                        fmpz_mul_2exp(e, fmpr_manref(fmprb_midref(y)), *exp_exp);
                        fmprb_pow_fmpz_binexp(z, x, e, prec);
                    }

                    fmpz_clear(e);
                    return;
                }
            }

            /* (-x)^n = (-1)^n * x^n */
            if (fmpz_sgn(exp_exp) >= 0 && fmpr_sgn(fmprb_midref(x)) < 0)
            {
                fmprb_t t;
                int odd;
                fmprb_init(t);
                fmprb_neg(t, x);
                odd = fmpz_is_zero(exp_exp);
                _fmprb_pow_exp(z, t, y, prec);
                if (odd)
                    fmprb_neg(z, z);
                fmprb_clear(t);
                return;
            }
        }
    }

    _fmprb_pow_exp(z, x, y, prec);
}

