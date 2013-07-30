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
#include "bernoulli.h"

double log2bern_approx(double n)
{
    return 1 + ((n+0.5)*log(n) - n - (n-0.5)*log(2*3.14159265358979323)) * (1. / log(2));
}

int main()
{
    long i, bound;
    double a, b;
    fmpq_t q;
    fmpr_t t;

    printf("bound_2exp_si....");
    fflush(stdout);

    fmpq_init(q);
    fmpr_init(t);

    for (i = 0; i < 1000; i++)
    {
        arith_bernoulli_number(q, i);
        bound = bernoulli_bound_2exp_si(i);

        fmpr_set_round_fmpz(t, fmpq_numref(q), 32, FMPR_RND_UP);
        fmpr_div_fmpz(t, t, fmpq_denref(q), 32, FMPR_RND_UP);

        if (fmpr_cmpabs_2exp_si(t, bound) > 0)
        {
            printf("FAIL: %ld\n", i);
            fmpr_print(t); printf("\n\n");
            printf("%ld\n", bound); printf("\n\n");
            abort();
        }
    }

    fmpq_clear(q);
    fmpr_clear(t);

    for (i = 100; i < 4000000; i += 1)
    {
        i += (i & 1);
        a = bernoulli_bound_2exp_si(i);
        b = log2bern_approx(i);

        if (a < b || a > 1.01 * b)
        {
            printf("FAIL: %ld\n", i);
            printf("%ld: %f %f %f\n", i, a, b, (float) a / b);
            abort();
        }
    }

    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

