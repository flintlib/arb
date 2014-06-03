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

    printf("d_log_lower_bound....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000; iter++)
    {
        mpfr_t t;
        double x, y, z;

        mpfr_init2(t, 53);

        x = d_randtest2(state);
        x = ldexp(x, 100 - n_randint(state, 200));

        switch (n_randint(state, 10))
        {
            case 0:
                x = 1.0 + x;
                break;
            case 1:
                x = 1.0 - x;
                break;
            case 2:
                x = D_INF;
                break;
            case 3:
                x = 0.0;
                break;
            case 4:
                x = D_NAN;
                break;
            default:
                break;
        }

        y = mag_d_log_lower_bound(x);

        mpfr_set_d(t, x, MPFR_RNDD);
        mpfr_log(t, t, MPFR_RNDD);
        z = mpfr_get_d(t, MPFR_RNDD);

        if (y > z || fabs(y-z) > 0.000001 * fabs(z))
        {
            printf("FAIL\n");
            printf("x = %.20g\n", x);
            printf("y = %.20g\n", y);
            printf("z = %.20g\n", z);
            abort();
        }

        mpfr_clear(t);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

