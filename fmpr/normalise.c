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

void _fmpr_normalise(fmpz_t man, fmpz_t exp, long prec, fmpr_rnd_t rnd)
{
    if (fmpz_is_zero(man))
    {
        /* this should perhaps crash, to avoid ambiguity? */
        fmpz_zero(exp);
        return;
    }
    else
    {
        long sign, bc, val, exp_shift;

        bc = fmpz_bits(man);
        sign = fmpz_sgn(man);
        exp_shift = 0;

        if (bc > prec)
        {
            exp_shift = bc - prec;

            if (rnd == FMPR_RND_NEAR)
            {
                abort();
            }
            else if (rnd == FMPR_RND_DOWN)
            {
                fmpz_tdiv_q_2exp(man, man, exp_shift);
            }
            else if (rnd == FMPR_RND_FLOOR)
            {
                fmpz_fdiv_q_2exp(man, man, exp_shift);
            }
            else if (rnd == FMPR_RND_CEIL)
            {
                fmpz_cdiv_q_2exp(man, man, exp_shift);
            }
            else
            {
                if (sign > 0)
                    fmpz_cdiv_q_2exp(man, man, exp_shift);
                else
                    fmpz_fdiv_q_2exp(man, man, exp_shift);
            }
        }

        val = fmpz_val2(man);
        exp_shift += val;

        if (val != 0)
            fmpz_tdiv_q_2exp(man, man, val);

        if (exp_shift != 0)
            fmpz_add_ui(exp, exp, exp_shift);
    }
}
