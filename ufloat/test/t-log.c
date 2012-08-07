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

#include <mpfr.h>
#include "ufloat.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("log....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000000; iter++)
    {
        mpfr_t x, y, z;
        ufloat_t u, v;

        if (n_randint(state, 2))
            ufloat_randtest(u, state, 10);
        else
            ufloat_randtest(u, state, 1L << 30);

        mpfr_init2(x, FLINT_BITS);
        mpfr_init2(y, FLINT_BITS);
        mpfr_init2(z, FLINT_BITS);

        ufloat_log(v, u);

        ufloat_get_mpfr(x, u);
        ufloat_get_mpfr(y, v);

        mpfr_log(z, x, MPFR_RNDU);

        if ((FLINT_BIT_COUNT(v->man) != UFLOAT_PREC) && (v->man != 0))
        {
            printf("wrong number of bits!\n");
            abort();
        }

        if (mpfr_cmp(z, y) > 0)
        {
            printf("fail!\n");
            printf("x = "); mpfr_printf("%.20Rg", x); printf("\n\n");
            printf("y = "); mpfr_printf("%.20Rg", y); printf("\n\n");
            printf("z = "); mpfr_printf("%.20Rg", z); printf("\n\n");
            abort();
        }

        mpfr_clear(x);
        mpfr_clear(y);
        mpfr_clear(z);
    }

    flint_randclear(state);
    mpfr_free_cache();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

