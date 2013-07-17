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

#include "zeta.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("series....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        fmpcb_t s, a;
        fmpcb_ptr z1, z2;
        long i, len, prec1, prec2;
        int deflate;

        fmpcb_init(s);
        fmpcb_init(a);

        if (n_randint(state, 2))
        {
            fmpcb_randtest(s, state, 1 + n_randint(state, 300), 3);
        }
        else
        {
            fmprb_set_ui(fmpcb_realref(s), 1);
            fmprb_mul_2exp_si(fmpcb_realref(s), fmpcb_realref(s), -1);
            fmprb_randtest(fmpcb_imagref(s), state, 1 + n_randint(state, 300), 4);
        }

        switch (n_randint(state, 3))
        {
            case 0:
                fmpcb_randtest(a, state, 1 + n_randint(state, 300), 3);
                break;
            case 1:
                fmprb_randtest(fmpcb_realref(a), state, 1 + n_randint(state, 300), 3);
                break;
            case 2:
                fmpcb_one(a);
                break;
        }

        prec1 = 2 + n_randint(state, 300);
        prec2 = prec1 + 30;
        len = 1 + n_randint(state, 20);

        deflate = n_randint(state, 2);

        z1 = _fmpcb_vec_init(len);
        z2 = _fmpcb_vec_init(len);

        zeta_series(z1, s, a, deflate, len, prec1);
        zeta_series(z2, s, a, deflate, len, prec2);

        for (i = 0; i < len; i++)
        {
            if (!fmpcb_overlaps(z1 + i, z2 + i))
            {
                printf("FAIL: overlap\n\n");
                printf("iter = %ld\n", iter);
                printf("deflate = %d, len = %ld, i = %ld\n\n", deflate, len, i);
                printf("s = "); fmpcb_printd(s, prec1 / 3.33); printf("\n\n");
                printf("a = "); fmpcb_printd(a, prec1 / 3.33); printf("\n\n");
                printf("z1 = "); fmpcb_printd(z1 + i, prec1 / 3.33); printf("\n\n");
                printf("z2 = "); fmpcb_printd(z2 + i, prec2 / 3.33); printf("\n\n");
                abort();
            }
        }

        fmpcb_clear(a);
        fmpcb_clear(s);
        _fmpcb_vec_clear(z1, len);
        _fmpcb_vec_clear(z2, len);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
