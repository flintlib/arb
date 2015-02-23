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

#include "arb.h"
#include "bernoulli.h"

void
arb_bernoulli_ui(arb_t b, ulong n, long prec)
{
    if (n < bernoulli_cache_num)
    {
        arb_set_fmpq(b, bernoulli_cache + n, prec);
    }
    else
    {
        int use_frac;

        use_frac = (n < BERNOULLI_SMALL_NUMER_LIMIT) || (n % 2 != 0);
        if (!use_frac && n < ULONG_MAX / 1000)
            use_frac = (prec > bernoulli_global_prec(n));

        if (use_frac)
        {
            fmpq_t t;
            fmpq_init(t);
            bernoulli_fmpq_ui(t, n);
            arb_set_fmpq(b, t, prec);
            fmpq_clear(t);
        }
        else
        {
            arb_bernoulli_ui_zeta(b, n, prec);
        }
    }
}

