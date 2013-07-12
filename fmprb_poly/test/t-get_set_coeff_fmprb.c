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

    Copyright (C) 2009 William Hart
    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "fmprb_poly.h"

int
main(void)
{
    int i, j, result;
    flint_rand_t state;

    printf("get/set_coeff_fmprb....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 100; i++)
    {
        fmprb_poly_t a;
        fmprb_t x1, x2;
        slong coeff, len;

        fmprb_poly_init(a);
        fmprb_init(x1);
        fmprb_init(x2);
        len = n_randint(state, 100) + 1;

        for (j = 0; j < 100; j++)
        {
            fmprb_randtest(x1, state, 2 + n_randint(state, 200), 10);
            coeff = n_randint(state, len);
            fmprb_poly_set_coeff_fmprb(a, coeff, x1);
            fmprb_poly_get_coeff_fmprb(x2, a, coeff);

            result = (fmprb_equal(x1, x2));
            if (!result)
            {
                printf("FAIL:\n");
                printf("x1 = "), fmprb_print(x1), printf("\n");
                printf("x2 = "), fmprb_print(x2), printf("\n");
                printf("coeff = %ld, length = %ld\n", coeff, len);
                abort();
            }
        }

        fmprb_clear(x1);
        fmprb_clear(x2);
        fmprb_poly_clear(a);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}

