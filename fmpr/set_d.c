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
fmpr_set_d(fmpr_t x, double v)
{
#if FLINT_BITS == 64
    mp_limb_t h, sign, exp, frac;
    long real_exp, real_man;
    union { double uf; mp_limb_t ul; } u;

    u.uf = v;
    h = u.ul;

    sign = h >> 63;
    exp = (h << 1) >> 53;
    frac = (h << 12) >> 12;

    if (exp == 0 && frac == 0)
    {
        fmpr_zero(x);
    }
    else if (exp == 0x7ff)
    {
        if (frac == 0)
        {
            if (sign)
                fmpr_neg_inf(x);
            else
                fmpr_pos_inf(x);
        }
        else
        {
            fmpr_nan(x);
        }
    }
    else
    {
        real_exp = exp - 1023 - 52;

        frac |= (1UL << 52);
        real_man = sign ? (-frac) : frac;

        fmpr_set_si_2exp_si(x, real_man, real_exp);
    }
#else
    mpfr_t t;
    mp_limb_t tmp[2];

    t->_mpfr_prec = 53;
    t->_mpfr_sign = 1;
    t->_mpfr_exp = 0;
    t->_mpfr_d = tmp;

    mpfr_set_d(t, v, MPFR_RNDD);

    fmpr_set_mpfr(x, t);
#endif
}

