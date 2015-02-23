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

    Copyright (C) 2015 Fredrik Johansson

******************************************************************************/

#include <string.h>
#include "arb.h"

int main()
{
    flint_rand_t state;

    printf("digits_round_inplace....");
    fflush(stdout);
    flint_randinit(state);

    {
        char s[30];
        long i, j, len, n;
        mp_bitcnt_t shift;
        fmpz_t inp, out, err, t;
        arf_rnd_t rnd;

        fmpz_init(inp);
        fmpz_init(out);
        fmpz_init(err);
        fmpz_init(t);

        for (i = 0; i < 100000; i++)
        {
            len = 1 + n_randint(state, 20);
            n = 1 + n_randint(state, 20);

            s[0] = (n_randint(state, 9) + '1');

            for (j = 1; j < len; j++)
                s[j] = (n_randint(state, 10) + '0');

            s[len] = '\0';

            fmpz_set_str(inp, s, 10);

            switch (n_randint(state, 3))
            {
                case 0:
                    rnd = ARF_RND_DOWN;
                    break;
                case 1:
                    rnd = ARF_RND_UP;
                    break;
                default:
                    rnd = ARF_RND_NEAR;
                    break;
            }

            _arb_digits_round_inplace(s, &shift, err, n, rnd);

            fmpz_set_str(out, s, 10);
            fmpz_set_ui(t, 10);
            fmpz_pow_ui(t, t, shift);
            fmpz_mul(t, t, out);
            fmpz_add(t, t, err);

            if (!fmpz_equal(t, inp) || (rnd == ARF_RND_UP && fmpz_sgn(err) > 0))
            {
                printf("FAIL!\n");
                printf("inp = "); fmpz_print(inp); printf("\n\n");
                printf("shift = %ld\n\n", shift);
                printf("err = "); fmpz_print(err); printf("\n\n");
                printf("out = "); fmpz_print(out); printf("\n\n");
                printf(" t  = "); fmpz_print(t); printf("\n\n");
                abort();
            }
        }

        fmpz_clear(inp);
        fmpz_clear(out);
        fmpz_clear(err);
        fmpz_clear(t);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

