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

#include <stdio.h>
#include <stdlib.h>
#include "bernoulli.h"
#include "ulong_extras.h"
#include "nmod_poly.h"
#include "nmod_vec.h"

int main()
{
    flint_rand_t state;
    long nmax, n, bound, count;
    mp_limb_t p, pinv, m1, m2;
    nmod_poly_t A;

    printf("rev....");
    fflush(stdout);
    flint_randinit(state);

    bound = 100000;

    p = n_nextprime(1UL << (FLINT_BITS - 1), 0);
    pinv = n_preinvert_limb(p);

    nmod_poly_init(A, p);
    nmod_poly_set_coeff_ui(A, 1, 1);
    nmod_poly_exp_series(A, A, bound);
    nmod_poly_shift_right(A, A, 1);
    nmod_poly_inv_series(A, A, bound);

    m1 = 1;
    for (n = 0; n < A->length; n++)
    {
        A->coeffs[n] = n_mulmod2_preinv(A->coeffs[n], m1, p, pinv);
        m1 = n_mulmod2_preinv(m1, n + 1, p, pinv);
    }

    for (nmax = 0; nmax < bound; nmax = 1.5 * nmax + 2)
    {
        fmpz_t numer, denom;
        bernoulli_rev_t iter;

        fmpz_init(numer);
        fmpz_init(denom);

        nmax += (nmax % 2);

        bernoulli_rev_init(iter, nmax);

        if (nmax < 8000)
            count = 4000;
        else
            count = 100;

        /* printf("nmax = %ld, count = %ld\n", nmax, count); */

        for (n = nmax; n >= 0 && count > 0; n -= 2, count--)
        {
            bernoulli_rev_next(numer, denom, iter);   

            m1 = fmpz_fdiv_ui(numer, p);
            m2 = fmpz_fdiv_ui(denom, p);
            m2 = n_invmod(m2, p);
            m1 = n_mulmod2_preinv(m1, m2, p, pinv);
            m2 = nmod_poly_get_coeff_ui(A, n);

            if (m1 != m2)
            {
                printf("FAIL:\n");
                printf("nmax = %ld, n = %ld\n", nmax, n);
                printf("m1 = %lu mod %lu\n", m1, p);
                printf("m2 = %lu mod %lu\n", m2, p);
                abort();
            }
        }

        bernoulli_rev_clear(iter);

        fmpz_clear(numer);
        fmpz_clear(denom);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

