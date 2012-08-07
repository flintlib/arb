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

#include "ufloat.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("set_ll_2exp....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        ufloat_t u;
        mp_limb_t hi, lo;
        fmpz_t a, b;
        long exp;

        hi = n_randtest(state);
        lo = n_randtest(state);
        exp = n_randint(state, 200);

        ufloat_set_ll_2exp(u, hi, lo, exp);

        fmpz_init(a);
        fmpz_init(b);

        fmpz_set_ui(a, hi);
        fmpz_mul_2exp(a, a, FLINT_BITS);
        fmpz_add_ui(a, a, lo);
        fmpz_mul_2exp(a, a, exp);

        ufloat_get_fmpz(b, u);

        if (fmpz_cmp(a, b) > 0 || fmpz_bits(b) > fmpz_bits(a) + 1)
        {
            printf("fail!\n");
            printf("hi = %lu\n", hi);
            printf("lo = %lu\n", lo);
            printf("exp = %ld\n", exp);
            printf("u = "); ufloat_print(u); printf("\n\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

