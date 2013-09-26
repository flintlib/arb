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

#include "fmprb_calc.h"

/* sin((pi/2)x) */
static int
sin_pi2_x(fmprb_ptr out, const fmprb_t inp, void * params, long order, long prec)
{
    fmprb_ptr x;

    x = _fmprb_vec_init(2);

    fmprb_set(x, inp);
    fmprb_one(x + 1);

    fmprb_const_pi(out, prec);
    fmprb_mul_2exp_si(out, out, -1);
    _fmprb_vec_scalar_mul(x, x, 2, out, prec);
    _fmprb_poly_sin_series(out, x, order, order, prec);

    _fmprb_vec_clear(x, 2);

    return 0;
}

int main()
{
    long iter;
    flint_rand_t state;

    printf("isolate_roots....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 40; iter++)
    {
        long m, r, a, b, maxdepth, maxeval, maxfound, prec, i, j, num;
        fmprb_ptr blocks;
        int * info;
        fmprb_t t, interval;
        fmpz_t nn;

        prec = 2 + n_randint(state, 50);

        m = n_randint(state, 80);
        r = 1 + n_randint(state, 80);
        a = m - r;
        b = m + r;

        maxdepth = 1 + n_randint(state, 60);
        maxeval = 1 + n_randint(state, 5000);
        maxfound = 1 + n_randint(state, 100);

        fmprb_init(interval);
        fmprb_init(t);
        fmpz_init(nn);

        fmpr_set_si(fmprb_midref(interval), m);
        fmpr_set_si(fmprb_radref(interval), r);

        num = fmprb_calc_isolate_roots(&blocks, &info, sin_pi2_x, NULL,
            interval, maxdepth, maxeval, maxfound, prec);

        /* check that all roots are accounted for */
        for (i = a; i <= b; i++)
        {
            if (i % 2 == 0)
            {
                int found = 0;

                for (j = 0; j < num; j++)
                {
                    if (fmprb_contains_si(blocks + j, i))
                    {
                        found = 1;
                        break;
                    }
                }

                if (!found)
                {
                    printf("FAIL: missing root %ld\n", i);
                    printf("a = %ld, b = %ld, maxdepth = %ld, maxeval = %ld, maxfound = %ld, prec = %ld\n",
                        a, b, maxdepth, maxeval, maxfound, prec);

                    for (j = 0; j < num; j++)
                    {
                        fmprb_printd(blocks + j, 15); printf("   %d \n", info[i]);
                    }

                    abort();
                }
            }
        }

        /* check that all reported single roots are good */
        for (i = 0; i < num; i++)
        {
            if (info[i] == 1)
            {
                /* b contains unique 2n -> b/2 contains unique n */
                fmprb_mul_2exp_si(t, blocks + i, -1);

                if (!fmprb_get_unique_fmpz(nn, t))
                {
                    printf("FAIL: bad root %ld\n", i);
                    printf("a = %ld, b = %ld, maxdepth = %ld, maxeval = %ld, maxfound = %ld, prec = %ld\n",
                        a, b, maxdepth, maxeval, maxfound, prec);

                    for (j = 0; j < num; j++)
                    {
                        fmprb_printd(blocks + j, 15); printf("   %d \n", info[i]);
                    }

                    abort();
                }
            }
        }

        _fmprb_vec_clear(blocks, num);
        flint_free(info);

        fmprb_clear(interval);
        fmprb_clear(t);
        fmpz_clear(nn);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

