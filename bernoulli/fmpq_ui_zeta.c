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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "bernoulli.h"
#include "arb.h"

void
_bernoulli_fmpq_ui_zeta(fmpz_t num, fmpz_t den, ulong n)
{
    slong prec;
    arb_t t;

    arith_bernoulli_number_denom(den, n);

    if (n % 2)
    {
        fmpz_set_si(num, -(n == 1));
        return;
    }

    if (n < BERNOULLI_SMALL_NUMER_LIMIT)
    {
        fmpz_set_si(num, _bernoulli_numer_small[n / 2]);
        return;
    }

    arb_init(t);

    for (prec = arith_bernoulli_number_size(n) + fmpz_bits(den) + 2; ; prec += 20)
    {
        arb_bernoulli_ui_zeta(t, n, prec);
        arb_mul_fmpz(t, t, den, prec);

        if (arb_get_unique_fmpz(num, t))
            break;

        printf("warning: %ld insufficient precision for Bernoulli number %lu\n", prec, n);
    }

    arb_clear(t);
}

