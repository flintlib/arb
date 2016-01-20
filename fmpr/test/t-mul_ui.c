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

    Copyright (C) 2012, 2014 Fredrik Johansson

******************************************************************************/

#include "fmpr.h"
#include "flint/ulong_extras.h"

static slong
fmpr_mul_ui_naive(fmpr_t z, const fmpr_t x, ulong y, slong prec, fmpr_rnd_t rnd)
{
    fmpr_t t; slong r;
    fmpr_init(t);
    fmpr_set_ui(t, y);
    r = fmpr_mul(z, x, t, prec, rnd);
    fmpr_clear(t);
    return r;
}

int main()
{
    slong iter, iter2;
    flint_rand_t state;

    flint_printf("mul_ui....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 3000; iter++)
    {
        fmpr_t x, z, v;
        ulong y;
        slong prec, r1, r2;
        fmpr_rnd_t rnd;

        fmpr_init(x);
        fmpr_init(z);
        fmpr_init(v);

        for (iter2 = 0; iter2 < 30; iter2++)
        {
            fmpr_randtest_special(x, state, 2000, 200);
            y = n_randtest(state);
            prec = 2 + n_randint(state, 2000);

            switch (n_randint(state, 4))
            {
                case 0:  rnd = FMPR_RND_DOWN; break;
                case 1:  rnd = FMPR_RND_UP; break;
                case 2:  rnd = FMPR_RND_FLOOR; break;
                default: rnd = FMPR_RND_CEIL; break;
            }

            switch (n_randint(state, 2))
            {
            case 0:
                r1 = fmpr_mul_ui(z, x, y, prec, rnd);
                r2 = fmpr_mul_ui_naive(v, x, y, prec, rnd);
                if (!fmpr_equal(z, v) || r1 != r2 || !fmpr_check_ulp(z, r1, prec))
                {
                    flint_printf("FAIL!\n");
                    flint_printf("x = "); fmpr_print(x); flint_printf("\n\n");
                    flint_printf("y = %wu\n\n", y);
                    flint_printf("z = "); fmpr_print(z); flint_printf("\n\n");
                    flint_printf("v = "); fmpr_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    abort();
                }
                break;

            default:
                fmpr_set(v, x);
                fmpr_set(z, x);
                r1 = fmpr_mul_ui(z, z, y, prec, rnd);
                r2 = fmpr_mul_ui_naive(v, v, y, prec, rnd);
                if (!fmpr_equal(z, v) || r1 != r2 || !fmpr_check_ulp(z, r1, prec))
                {
                    flint_printf("FAIL (aliasing 1)!\n");
                    flint_printf("x = "); fmpr_print(x); flint_printf("\n\n");
                    flint_printf("y = %wu\n\n", y);
                    flint_printf("z = "); fmpr_print(z); flint_printf("\n\n");
                    flint_printf("v = "); fmpr_print(v); flint_printf("\n\n");
                    flint_printf("r1 = %wd, r2 = %wd\n", r1, r2);
                    abort();
                }
                break;
            }
        }

        fmpr_clear(x);
        fmpr_clear(z);
        fmpr_clear(v);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
