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

    Copyright (C) 2015 Fredrik Johansson

******************************************************************************/

#include "acb_modular.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("elliptic_e....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 500; iter++)
    {
        acb_t m, w, K, Kp, E, Ep, r, pi2;
        long prec;

        acb_init(m);
        acb_init(w);
        acb_init(K);
        acb_init(Kp);
        acb_init(E);
        acb_init(Ep);
        acb_init(r);
        acb_init(pi2);

        prec = 2 + n_randint(state, 1000);
        acb_randtest(m, state, prec, 1 + n_randint(state, 100));

        acb_sub_ui(w, m, 1, prec);
        acb_neg(w, w);
        acb_const_pi(pi2, prec);
        acb_mul_2exp_si(pi2, pi2, -1);

        acb_modular_elliptic_k(K, m, prec);
        acb_modular_elliptic_k(Kp, w, prec);
        acb_modular_elliptic_e(E, m, prec);
        acb_modular_elliptic_e(Ep, w, prec);

        acb_mul(r, K, Ep, prec);
        acb_addmul(r, E, Kp, prec);
        acb_submul(r, K, Kp, prec);

        if (!acb_overlaps(r, pi2))
        {
            printf("FAIL (overlap)\n\n");

            printf("m = "); acb_printd(m, 30); printf("\n\n");
            printf("w = "); acb_printd(w, 30); printf("\n\n");
            printf("K = "); acb_printd(K, 30); printf("\n\n");
            printf("Kp = "); acb_printd(Kp, 30); printf("\n\n");
            printf("E = "); acb_printd(E, 30); printf("\n\n");
            printf("Ep = "); acb_printd(Ep, 30); printf("\n\n");
            printf("r = "); acb_printd(r, 30); printf("\n\n");
            abort();
        }

        acb_clear(m);
        acb_clear(w);
        acb_clear(K);
        acb_clear(Kp);
        acb_clear(E);
        acb_clear(Ep);
        acb_clear(r);
        acb_clear(pi2);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

