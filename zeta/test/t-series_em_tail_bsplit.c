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

#include "zeta.h"

void _fmpcb_poly_fmpcb_invpow_cpx(fmpcb_ptr res, const fmpcb_t N, const fmpcb_t c, long trunc, long prec);

int main()
{
    long iter;
    flint_rand_t state;

    printf("series_em_tail_bsplit....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        fmpcb_t Na, s;
        fmpcb_ptr z1, z2, Nasx;
        long i, M, len, prec;

        prec = 2 + n_randint(state, 400);
        len = 1 + n_randint(state, 30);
        M = n_randint(state, 40);

        fmpcb_init(Na);
        fmpcb_init(s);

        Nasx = _fmpcb_vec_init(len);
        z1 = _fmpcb_vec_init(len);
        z2 = _fmpcb_vec_init(len);

        fmpcb_randtest(Na, state, prec, 4);
        fmpcb_randtest(s, state, prec, 4);

        _fmpcb_poly_fmpcb_invpow_cpx(Nasx, Na, s, len, prec);

        zeta_em_tail_naive(z1, s, Na, Nasx, M, len, prec);
        zeta_em_tail_bsplit(z2, s, Na, Nasx, M, len, prec);

        for (i = 0; i < len; i++)
        {
            if (!fmpcb_overlaps(z1 + i, z2 + i))
            {
                printf("FAIL: overlap\n\n");
                printf("iter = %ld\n", iter);
                printf("prec = %ld, len = %ld, M = %ld\n", prec, len, M);
                printf("s = "); fmpcb_printd(s, prec / 3.33); printf("\n\n");
                printf("Na = "); fmpcb_printd(Na, prec / 3.33); printf("\n\n");
                printf("z1 = "); fmpcb_printd(z1 + i, prec / 3.33); printf("\n\n");
                printf("z2 = "); fmpcb_printd(z2 + i, prec / 3.33); printf("\n\n");
                abort();
            }
        }

        fmpcb_clear(Na);
        fmpcb_clear(s);
        _fmpcb_vec_clear(z1, len);
        _fmpcb_vec_clear(z2, len);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

