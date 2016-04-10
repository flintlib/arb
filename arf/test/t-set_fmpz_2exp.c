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

#include "arf.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("set_fmpz_2exp....");
    fflush(stdout);

    flint_randinit(state);

    /* test exact roundtrip R -> Q -> R */
    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        slong bits;
        arf_t x, z;
        fmpz_t y, e;

        bits = 2 + n_randint(state, 200);

        arf_init(x);
        arf_init(z);
        fmpz_init(y);
        fmpz_init(e);

        arf_randtest(x, state, bits, 1 + n_randint(state, 100));
        arf_randtest(z, state, bits, 1 + n_randint(state, 100));

        arf_get_fmpz_2exp(y, e, x);
        arf_set_fmpz_2exp(z, y, e);

        if (!arf_equal(x, z))
        {
            flint_printf("FAIL\n\n");
            flint_printf("bits: %wd\n", bits);
            flint_printf("x = "); arf_print(x); flint_printf("\n\n");
            flint_printf("y = "); fmpz_print(y); flint_printf("\n\n");
            flint_printf("e = "); fmpz_print(e); flint_printf("\n\n");
            flint_printf("z = "); arf_print(z); flint_printf("\n\n");
            abort();
        }

        arf_clear(x);
        arf_clear(z);
        fmpz_clear(y);
        fmpz_clear(e);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
