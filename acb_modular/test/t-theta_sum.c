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

#include "acb_modular.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("theta_sum....");
    fflush(stdout);

    flint_randinit(state);

    /* Very weak test, just testing the error bounds and not
       that we compute the right functions */
    for (iter = 0; iter < 10000; iter++)
    {
        acb_ptr t1a, t1b, t2a, t2b, t3a, t3b, t4a, t4b;
        acb_t w, q;
        int w_is_unit;
        long prec0, e0, prec1, prec2, len1, len2, i;

        acb_init(w);
        acb_init(q);

        e0 = 1 + n_randint(state, 100);
        prec0 = 2 + n_randint(state, 3000);
        prec1 = 2 + n_randint(state, 3000);
        prec2 = 2 + n_randint(state, 3000);
        len1 = 1 + n_randint(state, 30);
        len2 = 1 + n_randint(state, 30);

        t1a = _acb_vec_init(len1);
        t2a = _acb_vec_init(len1);
        t3a = _acb_vec_init(len1);
        t4a = _acb_vec_init(len1);

        t1b = _acb_vec_init(len2);
        t2b = _acb_vec_init(len2);
        t3b = _acb_vec_init(len2);
        t4b = _acb_vec_init(len2);

        if (n_randint(state, 2))
        {
            arb_randtest(acb_realref(q), state, prec0, e0);
            arb_zero(acb_imagref(q));
            acb_exp_pi_i(w, q, prec0);
            w_is_unit = n_randint(state, 2);
        }
        else
        {
            acb_randtest(w, state, prec0, e0);
            w_is_unit = 0;
        }

        acb_randtest(q, state, prec0, e0);

        for (i = 0; i < len1; i++)
        {
            acb_randtest(t1a + i, state, prec0, e0);
            acb_randtest(t2a + i, state, prec0, e0);
            acb_randtest(t3a + i, state, prec0, e0);
            acb_randtest(t4a + i, state, prec0, e0);
        }

        for (i = 0; i < len2; i++)
        {
            acb_randtest(t1b + i, state, prec0, e0);
            acb_randtest(t2b + i, state, prec0, e0);
            acb_randtest(t3b + i, state, prec0, e0);
            acb_randtest(t4b + i, state, prec0, e0);
        }

        acb_modular_theta_sum(t1a, t2a, t3a, t4a,
            w, w_is_unit, q, len1, prec1);

        acb_modular_theta_sum(t1b, t2b, t3b, t4b,
            w, w_is_unit & n_randint(state, 2), q, len2, prec2);

        for (i = 0; i < FLINT_MIN(len1, len2); i++)
        {
            if (!acb_overlaps(t1a + i, t1b + i)
                || !acb_overlaps(t2a + i, t2b + i)
                || !acb_overlaps(t3a + i, t3b + i)
                || !acb_overlaps(t4a + i, t4b + i))
            {
                printf("FAIL (overlap)  iter = %ld\n", iter);
                printf("len1 = %ld, len2 = %ld, prec1 = %ld, prec2 = %ld\n\n",
                    len1, len2, prec1, prec2);
                printf("i = %ld\n\n", i);
                printf("q = "); acb_printd(q, 50); printf("\n\n");
                printf("w = "); acb_printd(w, 50); printf("\n\n");
                printf("t1a = "); acb_printd(t1a + i, 50); printf("\n\n");
                printf("t1b = "); acb_printd(t1b + i, 50); printf("\n\n");
                printf("t2a = "); acb_printd(t2a + i, 50); printf("\n\n");
                printf("t2b = "); acb_printd(t2b + i, 50); printf("\n\n");
                printf("t3a = "); acb_printd(t3a + i, 50); printf("\n\n");
                printf("t3b = "); acb_printd(t3b + i, 50); printf("\n\n");
                printf("t4a = "); acb_printd(t4a + i, 50); printf("\n\n");
                printf("t4b = "); acb_printd(t4b + i, 50); printf("\n\n");
                abort();
            }
        }

        _acb_vec_clear(t1a, len1);
        _acb_vec_clear(t2a, len1);
        _acb_vec_clear(t3a, len1);
        _acb_vec_clear(t4a, len1);
        _acb_vec_clear(t1b, len2);
        _acb_vec_clear(t2b, len2);
        _acb_vec_clear(t3b, len2);
        _acb_vec_clear(t4b, len2);

        acb_clear(w);
        acb_clear(q);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

