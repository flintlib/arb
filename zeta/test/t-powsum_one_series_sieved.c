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

    printf("powsum_one_series_sieved....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 2000; iter++)
    {
        fmpcb_t s, a;
        fmpcb_ptr z1, z2;
        long i, n, len, prec;

        fmpcb_init(s);
        fmpcb_init(a);

        if (n_randint(state, 2))
        {
            fmpcb_randtest(s, state, 1 + n_randint(state, 200), 3);
        }
        else
        {
            fmprb_set_ui(fmpcb_realref(s), 1);
            fmprb_mul_2exp_si(fmpcb_realref(s), fmpcb_realref(s), -1);
            fmprb_randtest(fmpcb_imagref(s), state, 1 + n_randint(state, 200), 4);
        }

        fmpcb_one(a);

        prec = 2 + n_randint(state, 200);
        n = n_randint(state, 100);
        len = 1 + n_randint(state, 10);

        z1 = _fmpcb_vec_init(len);
        z2 = _fmpcb_vec_init(len);

        zeta_powsum_series_naive(z1, s, a, n, len, prec);
        zeta_powsum_one_series_sieved(z2, s, n, len, prec);

        for (i = 0; i < len; i++)
        {
            if (!fmpcb_overlaps(z1 + i, z2 + i))
            {
                printf("FAIL: overlap\n\n");
                printf("iter = %ld\n", iter);
                printf("n = %ld, prec = %ld, len = %ld, i = %ld\n\n", n, prec, len, i);
                printf("s = "); fmpcb_printd(s, prec / 3.33); printf("\n\n");
                printf("a = "); fmpcb_printd(a, prec / 3.33); printf("\n\n");
                printf("z1 = "); fmpcb_printd(z1 + i, prec / 3.33); printf("\n\n");
                printf("z2 = "); fmpcb_printd(z2 + i, prec / 3.33); printf("\n\n");
                abort();
            }
        }

        fmpcb_clear(a);
        fmpcb_clear(s);
        _fmpcb_vec_clear(z1, len);
        _fmpcb_vec_clear(z2, len);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
