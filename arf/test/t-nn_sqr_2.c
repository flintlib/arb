/*
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arf.h"

int main()
{
    slong ix, result;
    flint_rand_t state;

    flint_printf("nn_sqr_2....");
    fflush(stdout);

    flint_randinit(state);

    for (ix = 0; ix < 10000 * arb_test_multiplier(); ix++)
    {
        ulong r0, r1, r2, r3, s0, s1, s2, s3, a0, a1;

        a0 = n_randtest(state);
        a1 = n_randtest(state);

        nn_mul_2x2(r3, r2, r1, r0, a1, a0, a1, a0);
        nn_sqr_2(s3, s2, s1, s0, a1, a0);

        result = (s3 == r3 && s2 == r2 && s1 == r1 && s0 == r0);
        if (!result)
        {
            flint_printf("FAIL!\n");
            flint_abort();
        }
    }

    flint_randclear(state);
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
