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

#include <stdio.h>
#include <stdlib.h>
#include "elefun.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("exp_fixed_precomp....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        long prec, prec2;

        fmpz_t y, yerr, exponent, x, xerr;
        fmpr_t ya, yb, xa, xb, t, u;

        prec = 2 + n_randint(state, 1000);
        prec2 = prec + 30;

        fmpz_init(x);
        fmpz_init(xerr);
        fmpz_init(y);
        fmpz_init(yerr);
        fmpz_init(exponent);

        fmpr_init(ya);
        fmpr_init(yb);
        fmpr_init(xa);
        fmpr_init(xb);
        fmpr_init(t);
        fmpr_init(u);

        fmpz_randtest(x, state, prec + 20);
        fmpz_randtest_unsigned(xerr, state, prec + 10);

        fmpz_randtest(y, state, 1500);
        fmpz_randtest(yerr, state, 1500);
        fmpz_randtest(exponent, state, 1500);

        elefun_exp_fixed_precomp(y, yerr, exponent, x, xerr, prec);

        fmpr_set_fmpz(t, x);
        fmpr_mul_2exp_si(t, t, -prec);
        fmpr_set_fmpz(u, xerr);
        fmpr_mul_2exp_si(u, u, -prec);
        fmpr_sub(xa, t, u, FMPR_PREC_EXACT, FMPR_RND_FLOOR);
        fmpr_add(xb, t, u, FMPR_PREC_EXACT, FMPR_RND_CEIL);

        fmpr_set_fmpz(t, y);
        fmpr_mul_2exp_si(t, t, -prec);
        fmpr_mul_2exp_fmpz(t, t, exponent);
        fmpr_set_fmpz(u, yerr);
        fmpr_mul_2exp_si(u, u, -prec);
        fmpr_mul_2exp_fmpz(u, u, exponent);
        fmpr_sub(ya, t, u, FMPR_PREC_EXACT, FMPR_RND_FLOOR);
        fmpr_add(yb, t, u, FMPR_PREC_EXACT, FMPR_RND_CEIL);

        fmpr_exp(t, xa, prec2, FMPR_RND_FLOOR);
        fmpr_exp(u, xb, prec2, FMPR_RND_CEIL);

        if (!(fmpr_cmp(ya, t) < 0 && fmpr_cmp(t, yb) < 0))
        {
            printf("FAIL (bounds)\n");
            printf("prec = %ld\n", prec);
            printf("xa = "); fmpr_printd(xa, 50); printf("\n");
            printf("xb = "); fmpr_printd(xb, 50); printf("\n");
            printf("ya = "); fmpr_printd(ya, 50); printf("\n");
            printf("yb = "); fmpr_printd(yb, 50); printf("\n");
            printf("exp(xa) = "); fmpr_printd(t, 50); printf("\n");
            printf("exp(xb) = "); fmpr_printd(u, 50); printf("\n");
            abort();
        }

        fmpz_clear(x);
        fmpz_clear(xerr);
        fmpz_clear(y);
        fmpz_clear(yerr);
        fmpz_clear(exponent);

        fmpr_clear(ya);
        fmpr_clear(yb);
        fmpr_clear(xa);
        fmpr_clear(xb);
        fmpr_clear(t);
        fmpr_clear(u);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

