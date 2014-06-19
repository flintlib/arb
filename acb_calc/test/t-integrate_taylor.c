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

#include "acb_calc.h"

/* sin(x) */
int
sin_x(acb_ptr out, const acb_t inp, void * params, long order, long prec)
{
    int xlen = FLINT_MIN(2, order);

    acb_set(out, inp);
    if (xlen > 1)
        acb_one(out + 1);

    _acb_poly_sin_series(out, out, xlen, order, prec);
    return 0;
}

int main()
{
    long iter;
    flint_rand_t state;

    printf("integrate_taylor....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 150; iter++)
    {
        acb_t ans, res, a, b;
        arf_t inr, outr;
        double t;
        long goal, prec;

        acb_init(ans);
        acb_init(res);
        acb_init(a);
        acb_init(b);
        arf_init(inr);
        arf_init(outr);

        goal = 2 + n_randint(state, 300);
        prec = 2 + n_randint(state, 300);

        acb_randtest(a, state, 1 + n_randint(state, 200), 2);
        acb_randtest(b, state, 1 + n_randint(state, 200), 2);

        acb_cos(ans, a, prec);
        acb_cos(res, b, prec);
        acb_sub(ans, ans, res, prec);

        t = (1 + n_randint(state, 20)) / 10.0;
        arf_set_d(inr, t);
        arf_set_d(outr, t + (1 + n_randint(state, 20)) / 5.0);

        acb_calc_integrate_taylor(res, sin_x, NULL,
            a, b, inr, outr, goal, prec);

        if (!acb_overlaps(res, ans))
        {
            printf("FAIL! (iter = %ld)\n", iter);
            printf("prec = %ld, goal = %ld\n", prec, goal);
            printf("inr = "); arf_printd(inr, 15); printf("\n");
            printf("outr = "); arf_printd(outr, 15); printf("\n");
            printf("a = "); acb_printd(a, 15); printf("\n");
            printf("b = "); acb_printd(b, 15); printf("\n");
            printf("res = "); acb_printd(res, 15); printf("\n\n");
            printf("ans = "); acb_printd(ans, 15); printf("\n\n");
            abort();
        }

        acb_clear(ans);
        acb_clear(res);
        acb_clear(a);
        acb_clear(b);
        arf_clear(inr);
        arf_clear(outr);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

