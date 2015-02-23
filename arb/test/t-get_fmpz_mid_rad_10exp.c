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

#include "arb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("get_fmpz_mid_rad_10exp....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 5000; iter++)
    {
        arb_t x, y, t;
        fmpz_t mid, rad, exp;
        long n, prec;

        arb_init(x);
        arb_init(y);
        arb_init(t);
        fmpz_init(mid);
        fmpz_init(rad);
        fmpz_init(exp);

        arb_randtest_special(x, state, 1 + n_randint(state, 500), 1 + n_randint(state, 500));
        n = 1 + n_randint(state, 500);
        prec = 2 + n_randint(state, 1500);

        arb_get_fmpz_mid_rad_10exp(mid, rad, exp, x, n);

        arf_set_fmpz(arb_midref(y), mid);
        mag_set_fmpz(arb_radref(y), rad);
        arb_set_ui(t, 10);
        arb_pow_fmpz(t, t, exp, prec);
        arb_mul(y, y, t, prec);

        if (arb_is_finite(x) && !arb_is_zero(x) &&
            fmpz_sizeinbase(mid, 10) < n && fmpz_sizeinbase(rad, 10) < n)
        {
            printf("FAIL (too few digits):\n\n");
            printf("x = "); arb_printd(x, 50); printf("\n\n");
            printf("y = "); arb_printd(y, 50); printf("\n\n");
            printf("mid = "); fmpz_print(mid); printf("\n\n");
            printf("rad = "); fmpz_print(rad); printf("\n\n");
            printf("exp = "); fmpz_print(exp); printf("\n\n");
            abort();
        }

        if (arb_is_finite(x) && !arb_contains(y, x))
        {
            printf("FAIL (containment):\n\n");
            printf("x = "); arb_printd(x, 50); printf("\n\n");
            printf("y = "); arb_printd(y, 50); printf("\n\n");
            printf("mid = "); fmpz_print(mid); printf("\n\n");
            printf("rad = "); fmpz_print(rad); printf("\n\n");
            printf("exp = "); fmpz_print(exp); printf("\n\n");
            abort();
        }

        arb_clear(x);
        arb_clear(y);
        arb_clear(t);
        fmpz_clear(mid);
        fmpz_clear(rad);
        fmpz_clear(exp);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
