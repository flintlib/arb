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

void
ufloat_set_mpfr(ufloat_t u, const mpfr_t x)
{
    if (mpfr_zero_p(x))
    {
        u->man = 0;
        u->exp = 0;
    }
    else
    {
        mpz_t t;
        long exp, bits;

        mpz_init(t);
        exp = mpfr_get_z_2exp(t, x);
        mpz_abs(t, t);

        bits = mpz_sizeinbase(t, 2);

        if (bits > UFLOAT_PREC)
        {
            mpz_cdiv_q_2exp(t, t, bits - UFLOAT_PREC);
            exp += bits - UFLOAT_PREC;
        }

        u->exp = exp + UFLOAT_PREC;
        u->man = mpz_get_ui(t);
        ufloat_normalise(u);

        mpz_clear(t);
    }
}
