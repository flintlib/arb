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

#include "arb.h"

void
arb_atan_arf_via_mpfr(arb_t z, const arf_t x, long prec)
{
    mpfr_t t, u;
    int exact;

    mpfr_init2(t, 2 + arf_bits(x));
    mpfr_init2(u, prec);

    mpfr_set_emin(MPFR_EMIN_MIN);
    mpfr_set_emax(MPFR_EMAX_MAX);

    arf_get_mpfr(t, x, MPFR_RNDD);
    exact = (mpfr_atan(u, t, MPFR_RNDD) == 0);

    arf_set_mpfr(arb_midref(z), u);
    if (!exact)
        arf_mag_set_ulp(arb_radref(z), arb_midref(z), prec);

    mpfr_clear(t);
    mpfr_clear(u);
}

int main()
{
    long iter;
    flint_rand_t state;

    printf("atan_arf_bb....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 5000; iter++)
    {
        arb_t x, y, z;
        long prec, prec2;

        arb_init(x);
        arb_init(y);
        arb_init(z);

        prec = 2 + n_randint(state, 8000);

        arb_randtest(x, state, 1 + n_randint(state, 8000), 3);
        mag_zero(arb_radref(x));

        if (n_randint(state, 2))
            arb_mul_2exp_si(x, x, 1 + n_randint(state, 40));
        else
            arb_mul_2exp_si(x, x, -n_randint(state, 1.5 * prec));

        if (!arf_is_special(arb_midref(x)))
            prec2 = prec + 100 + 2 * (-ARF_EXP(arb_midref(x)));
        else
            prec2 = prec + 100;

        arb_atan_arf_via_mpfr(y, arb_midref(x), prec2);
        arb_atan_arf_bb(z, arb_midref(x), prec);

        if (!arb_contains(z, y))
        {
            printf("FAIL: containment\n\n");
            printf("prec = %ld\n\n", prec);
            printf("x = "); arb_printd(x, 50); printf("\n\n");
            printf("y = "); arb_printd(y, 50); printf("\n\n");
            printf("z = "); arb_printd(z, 50); printf("\n\n");
            abort();
        }

        if (arb_rel_accuracy_bits(z) < prec - 2)
        {
            printf("FAIL: poor accuracy\n\n");
            printf("prec = %ld,  acc = %ld\n\n", prec, arb_rel_accuracy_bits(z));
            printf("x = "); arb_printd(x, 50); printf("\n\n");
            printf("y = "); arb_printd(y, 50); printf("\n\n");
            printf("z = "); arb_printd(z, 50); printf("\n\n");
            abort();
        }

        arb_atan_arf_bb(x, arb_midref(x), prec);

        if (!arb_overlaps(x, z))
        {
            printf("FAIL: aliasing\n\n");
            abort();
        }

        arb_clear(x);
        arb_clear(y);
        arb_clear(z);
    }

    flint_randclear(state);
    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

