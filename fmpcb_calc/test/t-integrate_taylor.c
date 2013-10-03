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

#include "fmpcb_calc.h"

/* sin(x) */
int
sin_x(fmpcb_ptr out, const fmpcb_t inp, void * params, long order, long prec)
{
    int xlen = FLINT_MIN(2, order);

    fmpcb_set(out, inp);
    if (xlen > 1)
        fmpcb_one(out + 1);

    _fmpcb_poly_sin_series(out, out, xlen, order, prec);
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
        fmpcb_t ans, res, a, b;
        fmpr_t inr, outr;
        double t;
        long goal, prec;

        fmpcb_init(ans);
        fmpcb_init(res);
        fmpcb_init(a);
        fmpcb_init(b);
        fmpr_init(inr);
        fmpr_init(outr);

        goal = 2 + n_randint(state, 300);
        prec = 2 + n_randint(state, 300);

        fmpcb_randtest(a, state, 1 + n_randint(state, 200), 2);
        fmpcb_randtest(b, state, 1 + n_randint(state, 200), 2);

        fmpcb_cos(ans, a, prec);
        fmpcb_cos(res, b, prec);
        fmpcb_sub(ans, ans, res, prec);

        t = (1 + n_randint(state, 20)) / 10.0;
        fmpr_set_d(inr, t);
        fmpr_set_d(outr, t + (1 + n_randint(state, 20)) / 5.0);

        fmpcb_calc_integrate_taylor(res, sin_x, NULL,
            a, b, inr, outr, goal, prec);

        if (!fmpcb_overlaps(res, ans))
        {
            printf("FAIL! (iter = %ld)\n", iter);
            printf("prec = %ld, goal = %ld\n", prec, goal);
            printf("inr = "); fmpr_printd(inr, 15); printf("\n");
            printf("outr = "); fmpr_printd(outr, 15); printf("\n");
            printf("a = "); fmpcb_printd(a, 15); printf("\n");
            printf("b = "); fmpcb_printd(b, 15); printf("\n");
            printf("res = "); fmpcb_printd(res, 15); printf("\n\n");
            printf("ans = "); fmpcb_printd(ans, 15); printf("\n\n");
            abort();
        }

        fmpcb_clear(ans);
        fmpcb_clear(res);
        fmpcb_clear(a);
        fmpcb_clear(b);
        fmpr_clear(inr);
        fmpr_clear(outr);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

