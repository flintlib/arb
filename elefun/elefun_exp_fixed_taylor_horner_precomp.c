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

#include "elefun.h"

#define NUM_INVERSE_FACTORIALS 256
#define INVERSE_FACTORIALS_PREC 1024

__thread fmpz inverse_factorials[NUM_INVERSE_FACTORIALS];
__thread int inverse_factorials_init = 0;

void
compute_inverse_factorials()
{
    int i;
    fmpz_t t;
    fmpz_init(t);
    fmpz_one(t);
    fmpz_mul_2exp(t, t, INVERSE_FACTORIALS_PREC);
    for (i = 0; i < NUM_INVERSE_FACTORIALS; i++)
    {
        fmpz_init(inverse_factorials + i);
        fmpz_fac_ui(inverse_factorials + i, i);
        fmpz_tdiv_q(inverse_factorials + i, t, inverse_factorials + i);
    }
    fmpz_clear(t);
    inverse_factorials_init = 1;
}

void
elefun_exp_fixed_taylor_horner_precomp(fmpz_t y, fmpz_t yerr, const fmpz_t x, long n, long prec)
{
    if (n == 0 || prec > INVERSE_FACTORIALS_PREC || n > NUM_INVERSE_FACTORIALS)
    {
        abort();
    }
    else if (n == 1)  /* 1 */
    {
        fmpz_one(y);
        fmpz_mul_2exp(y, y, prec);
        fmpz_zero(yerr);
    }
    else if (n == 2)  /* 1 + x */
    {
        fmpz_one(y);
        fmpz_mul_2exp(y, y, prec);
        fmpz_add(y, y, x);
        fmpz_zero(yerr);
    }
    else if (n == 3)  /* 1 + x + x^2 / 2 */
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_one(y);
        fmpz_mul_2exp(y, y, prec);
        fmpz_add(y, y, x);
        fmpz_mul_tdiv_q_2exp(t, x, x, prec + 1);
        fmpz_add(y, y, t);
        fmpz_one(yerr);
        fmpz_clear(t);
    }
    else
    {
        fmpz_t t;
        long i;

        if (!inverse_factorials_init)
            compute_inverse_factorials();

        fmpz_init(t);

        fmpz_tdiv_q_2exp(y, inverse_factorials + n - 1, INVERSE_FACTORIALS_PREC - prec);

        for (i = n - 2; i >= 0; i--)
        {
            fmpz_mul_tdiv_q_2exp(y, y, x, prec);
            fmpz_tdiv_q_2exp(t, inverse_factorials + i, INVERSE_FACTORIALS_PREC - prec);
            fmpz_add(y, y, t);
        }

        fmpz_set_ui(yerr, 2 * n - 1);

        fmpz_clear(t);
    }
}

