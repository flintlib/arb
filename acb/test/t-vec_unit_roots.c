/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("vec_unit_roots....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100; iter++)
    {
        acb_ptr vec;
        fmpq_t q;
        acb_t t;
        slong k;
        slong prec = 53;

        vec = _acb_vec_init(iter);
        acb_init(t);
        fmpq_init(q);

        _acb_vec_unit_roots(vec, iter, prec);

        for (k = 0; k < iter; k++)
        {
            fmpq_set_si(q, 2 * k, iter);
            arb_sin_cos_pi_fmpq(acb_imagref(t), acb_realref(t), q, prec);

            if (!acb_overlaps(vec + k, t))
            {
                flint_printf("FAIL: overlap\n\n");
                flint_printf("n = %wu  k = %wd\n\n", iter, k);
                flint_printf("vec = "); acb_printn(vec + k, 30, 0); flint_printf("\n\n");
                flint_printf("t = "); acb_printn(t, 30, 0); flint_printf("\n\n");
                flint_abort();
            }
        }

        _acb_vec_clear(vec, iter);
        acb_clear(t);
        fmpq_clear(q);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

