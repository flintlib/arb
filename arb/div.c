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


void _arb_div_shiftup(fmpz_t cmid, fmpz_t crad, fmpz_t cexp,
    const fmpz_t amid, const fmpz_t arad, const fmpz_t aexp,
    const fmpz_t bmid, const fmpz_t brad, const fmpz_t bexp,
    long shift)
{
    fmpz_t t, u, v;

    fmpz_init(t);
    fmpz_init(u);
    fmpz_init(v);

    /* a/b - (a+r)/(b-s) = r/(s-b) + a*s/(b*(s-b)) */

    /* t = minimize |s-b| */
    if (fmpz_sgn(bmid) >= 0)
    {
        fmpz_sub(t, bmid, brad);
    }
    else
    {
        fmpz_add(t, bmid, brad);
        fmpz_neg(t, t);
    }

    /* first terror term: r/|s-b| */
    fmpz_mul_2exp(u, arad, shift);
    fmpz_cdiv_q(u, u, t);

    /* second error term (likely small): a*s/(b*(s-b)) */
    fmpz_mul_2exp(v, amid, shift);
    fmpz_mul(v, v, brad);
    fmpz_abs(v, v);
    fmpz_cdiv_q_2exp(v, v, fmpz_bits(bmid) - 1);   /* fast approximate division */
    fmpz_cdiv_q_2exp(v, v, fmpz_bits(t) - 1);      /* fast approximate division */

    /* add error terms */
    fmpz_add(crad, u, v);

    /* compute quotient */
    fmpz_mul_2exp(v, amid, shift);
    fmpz_tdiv_q(cmid, v, bmid);

    /* error from main division */
    fmpz_add_ui(crad, crad, 1UL);

    /* subtract exponents */
    fmpz_sub(cexp, aexp, bexp);
    fmpz_sub_ui(cexp, cexp, shift);

    fmpz_clear(t);
    fmpz_clear(u);
    fmpz_clear(v);
}

void
arb_div(arb_t c, const arb_t a, const arb_t b)
{
    long abits, bbits, shift;

    if (arb_contains_zero(b))
    {
        printf("arb_div: division by zero\n");
        abort();
    }

    if (fmpz_is_zero(arb_radref(b)))
    {
        if (fmpz_is_zero(arb_radref(a)))
        {
            fmpz_t t;
            fmpz_init(t);
            fmpz_sub(t, arb_expref(a), arb_expref(b));
            _arb_set_fmpq(c, arb_midref(a), arb_midref(b));
            fmpz_add(arb_expref(c), arb_expref(c), t);
            fmpz_clear(t);
            return;
        }

        /* TODO: (a + r) / b = a / b + r / b */
    }


    abits = fmpz_bits(arb_midref(a));
    bbits = fmpz_bits(arb_midref(b));
    shift = arb_prec(c) - (abits - bbits);

    /* shift = FLINT_MAX(shift, 0);  more precise, but can be much slower */

    if (shift >= 0)
    {
        _arb_div_shiftup(arb_midref(c), arb_radref(c), arb_expref(c),
            arb_midref(a), arb_radref(a), arb_expref(a),
            arb_midref(b), arb_radref(b), arb_expref(b), shift);
    }
    else
    {
        fmpz_t amid, arad, aexp;

        fmpz_init(amid);
        fmpz_init(arad);
        fmpz_init(aexp);

        fmpz_add_ui(aexp, arb_expref(a), -shift);
        fmpz_tdiv_q_2exp(amid, arb_midref(a), -shift);
        fmpz_cdiv_q_2exp(arad, arb_radref(a), -shift);
        fmpz_add_ui(arad, arad, 1UL);

        _arb_div_shiftup(arb_midref(c), arb_radref(c), arb_expref(c),
            amid, arad, aexp, arb_midref(b), arb_radref(b), arb_expref(b), 0);

        fmpz_clear(amid);
        fmpz_clear(arad);
        fmpz_clear(aexp);
    }
}
