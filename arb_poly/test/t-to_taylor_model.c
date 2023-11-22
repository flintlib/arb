/*
    Copyright (C) 2022 Erik Postma

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("to_taylor_model....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 50 * arb_test_multiplier(); ++iter)
    {
        slong prec, glen, n_iterations, y_prec, i;
        arb_poly_t fp;
        arb_t x;
        arb_ptr g;
        mag_t radius;

        prec = 2 + n_randint(state, 500);

        arb_poly_init(fp);
        arb_init(x);
        mag_init(radius);

        arb_poly_randtest(fp, state, 1 + n_randint(state, 10), 1 + n_randint(state, 500), 10);
        arb_randtest(x, state, 1 + n_randint(state, 100), 10);

        glen = 1 + n_randint(state, 10); /* We can have glen >= len. That's not very useful, but it
                                          * should work. */
        g = _arb_vec_init(glen);

        _arb_poly_to_taylor_model(g, radius, fp->coeffs, fp->length, x, glen, prec);

        if(arb_is_exact(x))
        {
            n_iterations = 1;
            y_prec = prec + 30;
        }
        else
        {
            n_iterations = 25;
            y_prec = prec;
        }

        for(i = 0; i < n_iterations; ++i)
        {
            /* Generate a random point y = arb_midref(x) + yy inside the interval represented by
             * x. Verify that f(y) is contained in g(y) +/- radius. */
            arb_t y, yy, fy, gy;

            arb_init(y);
            arb_init(yy);
            arb_init(fy);
            arb_init(gy);

            if(arb_is_exact(x))
            {
                arb_zero(yy);
                arb_set(y, x);
            }
            else
            {
                arf_t rad_x_as_arf_t;
                arf_init(rad_x_as_arf_t);

                arb_urandom(yy, state, y_prec);

                arf_set_mag(rad_x_as_arf_t, arb_radref(x));
                arb_mul_arf(yy, yy, rad_x_as_arf_t, y_prec);
                arb_add_arf(y, yy, arb_midref(x), y_prec);

                arf_clear(rad_x_as_arf_t);
            }

            arb_poly_evaluate(fy, fp, y, y_prec + 30);

            _arb_poly_evaluate(gy, g, glen, yy, y_prec);
            flint_printf("gy before = "); arb_printd(gy, 15); flint_printf("\n");
            arb_add_error_mag(gy, radius);
            flint_printf("gy after  = "); arb_printd(gy, 15); flint_printf("\n");

            if(! (glen >= fp->length ? arb_overlaps(gy, fy) : arb_contains(gy, fy)))
            {
                flint_printf("FAIL\n\n");

                flint_printf("prec = %d\n", prec);
                flint_printf("len = %d\n", fp->length);
                flint_printf("glen = %d\n\n", glen);
                
                flint_printf("f = ");      arb_poly_printd(fp, 15); flint_printf("\n");
                flint_printf("x = ");      arb_printd(x, 15);       flint_printf("\n");
                flint_printf("g = [");
                for(i = 0; i < glen; ++i)
                {
                    flint_printf("(");
                    arb_printd(g + i, 15);
                    if(i < glen-1)
                        flint_printf("), ");
                    else
                        flint_printf(")");
                }
                flint_printf("]\n");
                flint_printf("radius = "); mag_printd(radius, 15);  flint_printf("\n\n");

                flint_printf("yy = ");     arb_printd(yy, 15);      flint_printf("\n");
                flint_printf("y = ");      arb_printd(y, 15);       flint_printf("\n");
                flint_printf("fy = ");     arb_printd(fy, 50);      flint_printf("\n");
                flint_printf("gy = ");     arb_printd(gy, 50);      flint_printf("\n\n");

                flint_abort();
            }

            arb_clear(gy);
            arb_clear(fy);
            arb_clear(y);
            arb_clear(yy);
        }
        
        mag_clear(radius);
        arb_clear(x);
        arb_poly_clear(fp);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
