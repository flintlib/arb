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

long
_fmpr_normalise_naive(fmpz_t man, fmpz_t exp, long prec, fmpr_rnd_t rnd)
{
    /* TODO: this should perhaps raise an exception to avoid ambiguity */
    if (fmpz_is_zero(man))
    {
        fmpz_zero(exp);
        return FMPR_RESULT_EXACT;
    }
    else
    {
        long bc, val;

        bc = fmpz_bits(man);
        val = fmpz_val2(man);

        if (bc - val <= prec)
        {
            if (val != 0)
            {
                fmpz_tdiv_q_2exp(man, man, val);
                fmpz_add_ui(exp, exp, val);
            }

            return FMPR_RESULT_EXACT;
        }
        else
        {
            long exp_shift = bc - prec;

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
                if (fmpz_sgn(man) > 0)
                    fmpz_cdiv_q_2exp(man, man, exp_shift);
                else
                    fmpz_fdiv_q_2exp(man, man, exp_shift);
            }

            val = fmpz_val2(man);
            exp_shift += val;

            if (val != 0)
                fmpz_tdiv_q_2exp(man, man, val);

            fmpz_add_ui(exp, exp, exp_shift);
            return val - (val == prec);
        }
    }
}
