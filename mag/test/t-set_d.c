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

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include "double_extras.h"
#include "mag.h"

/* XXX: d_randtest is not good enough */

#define EXP_MINUS_32 2.3283064365386962891e-10
#define EXP_MINUS_64 5.42101086242752217e-20

double
d_randtest2(flint_rand_t state)
{
    mp_limb_t m1, m2;
    double t;

    if (FLINT_BITS == 64)
    {
        m1 = n_randtest(state) | (UWORD(1) << (FLINT_BITS - 1));

        t = ((double) m1) * EXP_MINUS_64;
    }
    else
    {
        m1 = n_randtest(state) | (UWORD(1) << (FLINT_BITS - 1));
        m2 = n_randtest(state);

        t = ((double) m1) * EXP_MINUS_32 +
            ((double) m2) * EXP_MINUS_64;
    }

    return t;
}

int main()
{
    long iter;
    flint_rand_t state;

    printf("set_d....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        fmpr_t a, b, c;
        mag_t m;
        double x;

        fmpr_init(a);
        fmpr_init(b);
        fmpr_init(c);
        mag_init(m);

        x = d_randtest2(state);
        x = ldexp(x, 100 - n_randint(state, 200));

        if (n_randint(state, 100) == 0)
            x = 0.0;

        fmpr_set_d(a, x);
        mag_set_d(m, x);

        mag_get_fmpr(b, m);

        fmpr_set(c, a);
        fmpr_mul_ui(c, c, 1025, MAG_BITS, FMPR_RND_UP);
        fmpr_mul_2exp_si(c, c, -10);

        MAG_CHECK_BITS(m)

        if (!(fmpr_cmpabs(a, b) <= 0 && fmpr_cmpabs(b, c) <= 0))
        {
            printf("FAIL\n\n");
            printf("a = "); fmpr_print(a); printf("\n\n");
            printf("b = "); fmpr_print(b); printf("\n\n");
            printf("c = "); fmpr_print(c); printf("\n\n");
            abort();
        }

        fmpr_clear(a);
        fmpr_clear(b);
        fmpr_clear(c);
        mag_clear(m);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

