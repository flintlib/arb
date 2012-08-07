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

    printf("set_mpfr....");
    fflush(stdout);

    flint_randinit(state);

    /* Exact roundtrip */
    for (iter = 0; iter < 100000; iter++)
    {
        ufloat_t u, v;
        mpfr_t x;

        mpfr_init2(x, UFLOAT_PREC);

        ufloat_randtest(u, state, 10);

        ufloat_get_mpfr(x, u);
        ufloat_set_mpfr(v, x);

        if (!ufloat_equal(u, v))
        {
            printf("fail!\n");
            printf("u = "); ufloat_print(u); printf("\n\n");
            printf("v = "); ufloat_print(v); printf("\n\n");
            abort();
        }

        mpfr_clear(x);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

