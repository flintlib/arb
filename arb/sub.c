/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "arb.h"

void
arb_sub(arb_t c, const arb_t a, const arb_t b)
{
    fmpz e, f;

    e = *arb_expref(a);
    f = *arb_expref(b);

    /* both exponents equal and small (e.g. small integer case) */
    if (e == f)
    {
        /* subtract midpoints */
        fmpz_sub(arb_midref(c), arb_midref(a), arb_midref(b));

        /* add errors */
        if (fmpz_is_zero(arb_radref(a)) && fmpz_is_zero(arb_radref(b)))
            fmpz_zero(arb_radref(c));
        else
            fmpz_add(arb_radref(c), arb_radref(a), arb_radref(b));

        /* set exponent */
        _fmpz_set_si_small(arb_expref(c), e);
    }
    /* both exponents are small */
    else if (!COEFF_IS_MPZ(e) && !COEFF_IS_MPZ(f))
    {
        fmpz_t t;
        fmpz_init(t);

        /* TODO: be smart when one argument (or both) has rad = 0 */

        if (e >= f)
        {
            fmpz_tdiv_q_2exp(t, arb_midref(b), e - f);
            fmpz_sub(arb_midref(c), arb_midref(a), t);

            fmpz_cdiv_q_2exp(t, arb_radref(b), e - f);
            fmpz_add(arb_radref(c), arb_radref(a), t);
            fmpz_add_ui(arb_radref(c), arb_radref(c), 1UL);

            _fmpz_set_si_small(arb_expref(c), e);
        }
        else
        {
            fmpz_tdiv_q_2exp(t, arb_midref(a), f - e);
            fmpz_sub(arb_midref(c), t, arb_midref(b));

            fmpz_cdiv_q_2exp(t, arb_radref(a), f - e);
            fmpz_add(arb_radref(c), arb_radref(b), t);
            fmpz_add_ui(arb_radref(c), arb_radref(c), 1UL);

            _fmpz_set_si_small(arb_expref(c), f);
        }

        fmpz_clear(t);
    }
    else
    {
        /* not implemented */
        abort();
    }
}
