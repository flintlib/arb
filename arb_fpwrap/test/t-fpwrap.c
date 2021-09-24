/*
    Copyright (C) 2021 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "flint/double_extras.h"
#include "arb_fpwrap.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("fpwrap....");
    fflush(stdout);

    flint_randinit(state);

    /* correct rounding test */
    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        mpfr_t t;
        double x, y, z;

        mpfr_init2(t, 53);

        x = d_randtest(state) + n_randint(state, 100);

        mpfr_set_d(t, x, MPFR_RNDN);

        switch (n_randint(state, 4))
        {
            case 0:
                arb_fpwrap_double_log1p(&y, x, FPWRAP_CORRECT_ROUNDING);
                mpfr_log1p(t, t, MPFR_RNDN);
                break;
            case 1:
                arb_fpwrap_double_sqrt(&y, x, FPWRAP_CORRECT_ROUNDING);
                mpfr_sqrt(t, t, MPFR_RNDN);
                break;
            case 2:
                arb_fpwrap_double_exp(&y, x, FPWRAP_CORRECT_ROUNDING);
                mpfr_exp(t, t, MPFR_RNDN);
                break;
            default:
                arb_fpwrap_double_sin(&y, x, FPWRAP_CORRECT_ROUNDING);
                mpfr_sin(t, t, MPFR_RNDN);
                break;
        }

        z = mpfr_get_d(t, MPFR_RNDN);

        if (z != y)
        {
            flint_printf("FAIL: correct rounding\n\n");
            flint_abort();
        }

        mpfr_clear(t);
    }

    {
        double a[1], b[2];
        double z, y;

        z = 1.75;
        a[0] = 0.25;
        b[0] = 1.5;
        b[1] = -2.125;

        arb_fpwrap_double_hypgeom_pfq(&y, a, 1, b, 2, z, 0, 0);

        if (fabs(y - 0.68910385124070327187) > 1e-16)
        {
            flint_printf("FAIL: value 1\n\n");
            flint_abort();
        }

        arb_fpwrap_double_hypgeom_pfq(&y, a, 1, b, 2, z, 1, 0);

        if (fabs(y - (-0.21324224371323783595)) > 1e-16)
        {
            flint_printf("FAIL: value 2\n\n");
            flint_abort();
        }
    }

    {
        complex_double y, z;

        arb_fpwrap_cdouble_zeta_zero(&y, 1, FPWRAP_CORRECT_ROUNDING);
        arb_fpwrap_cdouble_zeta(&z, y, 0);

        if (fabs(z.real - (-1.0483650805588237388e-16)) > 1e-31)
        {
            flint_printf("FAIL: value 3\n\n");
            flint_abort();
        }

        if (fabs(z.imag - 6.5852592776051578103e-16) > 1e-31)
        {
            flint_printf("FAIL: value 4\n\n");
            flint_abort();
        }
    }

    {
        complex_double x, y;

        x.real = 1.0;
        x.imag = 1e-100;

        arb_fpwrap_cdouble_erf(&y, x, FPWRAP_ACCURATE_PARTS);

        if (fabs(y.imag - 4.1510749742059471164e-101) > 1e-116)
        {
            flint_printf("FAIL: value 5\n\n");
            flint_abort();
        }
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

