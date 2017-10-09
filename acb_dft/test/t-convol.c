/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dft.h"

typedef void (*do_f) (acb_ptr z, acb_srcptr x, acb_srcptr y, slong len, slong prec);

void
check_vec_eq_prec(acb_srcptr w1, acb_srcptr w2, slong len, slong prec, slong digits, ulong q, char f1[], char f2[])
{
    slong i;

    for (i = 0; i < len; i++)
    {
        if (!acb_overlaps(w1 + i, w2 + i))
        {
            flint_printf("FAIL\n\n");
            flint_printf("q = %wu, size = %wu\n", q, len);
            flint_printf("\nDFT differ from index %ld / %ld \n", i, len);
            flint_printf("\n%s =\n", f1);
            acb_vec_printd_index(w1, len, digits);
            flint_printf("\n%s =\n", f2);
            acb_vec_printd_index(w2, len, digits);
            flint_printf("\n\n");
            abort();
        }
        else if (!acb_is_zero(w1+i) && (acb_rel_accuracy_bits(w1 + i) < 30
                || acb_rel_accuracy_bits(w2 + i) < 30))
        {
            flint_printf("FAIL\n\n");
            flint_printf("q = %wu\n", q);
            flint_printf("\nDFT inaccurate from index %ld / %ld \n", i, len);
            flint_printf("\n%s =\n", f1);
            acb_printd(w1 + i, digits);
            flint_printf("\n%s =\n", f2);
            acb_printd(w2 + i, digits);
            flint_printf("\nerrors %ld & %ld [prec = %wu]\n",
                    acb_rel_accuracy_bits(w1 + i),
                    acb_rel_accuracy_bits(w2 + i), prec);
            abort();
        }
    }
}


int main()
{

    slong k;
    slong prec = 100, digits = 30;
    slong nq = 15;
    ulong q[15] = { 2, 3, 4, 5, 6, 23, 10, 15, 30, 59, 256, 308, 335, 344, 961};
    flint_rand_t state;

    slong f, nf = 3;
    do_f func[3] = { acb_dft_convol_naive, acb_dft_convol_rad2, acb_dft_convol_dft };
    char * name[3] = { "naive", "rad2", "dft" };

    flint_printf("convol....");
    fflush(stdout);

    flint_randinit(state);

    for (k = 0; k < nq; k++)
    {
        slong i;
        acb_ptr z1, z2, x, y;

        z1 = _acb_vec_init(q[k]);
        z2 = _acb_vec_init(q[k]);
        x = _acb_vec_init(q[k]);
        y = _acb_vec_init(q[k]);

        for (i = 0; i < q[k]; i++)
        {
            acb_set_si(x + i, q[k] - i);
            acb_set_si(y + i, i * i);
            /*
            acb_set_si(x + i, n_randint(state, q[k]));
            acb_set_si(y + i, n_randint(state, q[k]));
            */
        }

        for (f = 0; f < nf; f++)
        {

            acb_ptr z = (f == 0) ? z1 : z2;

            func[f](z, x, y, q[k], prec);

            if (f == 0)
                continue;

            check_vec_eq_prec(z1, z2, q[k], prec, digits, q[k], name[0], name[f]);

        }

        _acb_vec_clear(x, q[k]);
        _acb_vec_clear(y, q[k]);
        _acb_vec_clear(z1, q[k]);
        _acb_vec_clear(z2, q[k]);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}


