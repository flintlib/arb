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

static const double answers[10] = {
  1.04570093561423094, 2.0358667496686487, 4.82706400405656566,
  11.2347033850581819, 27.1331828778516522, 67.210305990439042,
  168.564351176369878, 427.496180295369909, 1093.56959348921777,
  2815.70144392142227
};

int main()
{
    long iter;
    flint_rand_t state;

    printf("cauchy_bound....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100; iter++)
    {
        fmprb_t b, radius, ans;
        fmpcb_t x;
        long r, prec, maxdepth;

        fmprb_init(b);
        fmprb_init(radius);
        fmprb_init(ans);
        fmpcb_init(x);

        fmpcb_set_ui(x, 5);

        r = 1 + n_randint(state, 10);
        fmprb_set_ui(radius, r);

        prec = 2 + n_randint(state, 100);
        maxdepth = n_randint(state, 10);

        fmpcb_calc_cauchy_bound(b, sin_x, NULL, x, radius, maxdepth, prec);

        fmpr_set_d(fmprb_midref(ans), answers[r-1]);
        fmpr_set_d(fmprb_radref(ans), 1e-8);

        if (!fmprb_overlaps(b, ans))
        {
            printf("FAIL\n");
            printf("r = %ld, prec = %ld, maxdepth = %ld\n\n", r, prec, maxdepth);
            fmprb_printd(b, 15); printf("\n\n");
            fmprb_printd(ans, 15); printf("\n\n");
            abort();
        }

        fmprb_clear(b);
        fmprb_clear(radius);
        fmprb_clear(ans);
        fmpcb_clear(x);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

