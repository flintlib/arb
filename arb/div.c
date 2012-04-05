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
arb_div(arb_t c, const arb_t a, const arb_t b)
{
    fmpz_t t, u, v;
    long abits, bbits, sshift, shift;

    if (arb_contains_zero(b))
    {
        printf("arb_div: division by zero\n");
        abort();
    }

    if (fmpz_is_zero(arb_radref(b)))
    {
        if (fmpz_is_zero(arb_radref(a)))
        {
            fmpz_init(t);
            fmpz_sub(t, arb_expref(a), arb_expref(b));
            _arb_set_fmpq(c, arb_midref(a), arb_midref(b));
            fmpz_add(arb_expref(c), arb_expref(c), t);
            fmpz_clear(t);
            return;
        }

        /* TODO: (a + r) / b = a / b + r / b */
    }

    fmpz_init(t);
    fmpz_init(u);
    fmpz_init(v);

    abits = fmpz_bits(arb_midref(a));
    bbits = fmpz_bits(arb_midref(b));
    sshift = arb_prec(c) - (abits - bbits);
    shift = FLINT_MAX(0, sshift);

    /* a/b - (a+r)/(b-s) = r/(s-b) + a*s/(b*(s-b)) */

    /* t = minimize |s-b| */
    if (fmpz_sgn(arb_midref(b)) >= 0)
    {
        fmpz_sub(t, arb_midref(b), arb_radref(b));
    }
    else
    {
        fmpz_add(t, arb_midref(b), arb_radref(b));
        fmpz_neg(t, t);
    }

    /* first terror term: r/|s-b| */
    fmpz_mul_2exp(u, arb_radref(a), shift);
    fmpz_cdiv_q(u, u, t);

    /* second error term (likely small): a*s/(b*(s-b)) */
    fmpz_mul_2exp(v, arb_midref(a), shift);
    fmpz_mul(v, v, arb_radref(b));
    fmpz_abs(v, v);
    fmpz_cdiv_q_2exp(v, v, bbits - 1);         /* fast approximate division */
    fmpz_cdiv_q_2exp(v, v, fmpz_bits(t) - 1);  /* fast approximate division */

    /* add error terms */
    fmpz_add(arb_radref(c), u, v);

    /* compute quotient; TODO: normalise afterwards if sshift < 0 */
    fmpz_mul_2exp(v, arb_midref(a), shift);
    fmpz_tdiv_q(arb_midref(c), v, arb_midref(b));

    /* error from main division */
    fmpz_add_ui(arb_radref(c), arb_radref(c), 1UL);

    /* subtract exponents */
    fmpz_sub(arb_expref(c), arb_expref(a), arb_expref(b));
    fmpz_sub_ui(arb_expref(c), arb_expref(c), shift);

    fmpz_clear(t);
    fmpz_clear(u);
    fmpz_clear(v);
}
