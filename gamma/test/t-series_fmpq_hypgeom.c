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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "fmprb.h"
#include "gamma.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("series_fmpq_hypgeom....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 250; iter++)
    {
        fmprb_ptr u, v;
        fmpq_t a;
        ulong p, q;
        long i, len, prec1, prec2;

        len = 1 + n_randint(state, 10);
        prec1 = 2 + n_randint(state, 1000);
        prec2 = prec1 + 30;

        u = _fmprb_vec_init(len);
        v = _fmprb_vec_init(len);
        fmpq_init(a);

        q = 1 + n_randint(state, 20);
        p = q - n_randint(state, q);
        fmpq_set_si(a, p, q);

        gamma_series_fmpq_hypgeom(u, a, len, prec1);
        gamma_series_fmpq_hypgeom(v, a, len, prec2);

        for (i = 0; i < len; i++)
        {
            if (!fmprb_overlaps(u + i, v + i))
            {
                printf("FAIL: overlap\n\n");
                printf("p = %lu, q = %lu, len = %ld, i = %ld\n\n", p, q, len, i);
                printf("u = "); fmprb_printd(u + i, prec1 / 3.33); printf("\n\n");
                printf("v = "); fmprb_printd(v + i, prec2 / 3.33); printf("\n\n");
                abort();
            }
        }

        fmpq_clear(a);
        _fmprb_vec_clear(u, len);
        _fmprb_vec_clear(v, len);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

