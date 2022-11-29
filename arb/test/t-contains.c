/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"

int
_fmpq_cmp2(const fmpz_t p, const fmpz_t q, const fmpz_t r, const fmpz_t s)
{
    int s1, s2, res;
    flint_bitcnt_t bp, bq, br, bs;
    fmpz_t t, u;

    if (!COEFF_IS_MPZ(*p) && !COEFF_IS_MPZ(*q) && !COEFF_IS_MPZ(*r) && !COEFF_IS_MPZ(*s))
    {
        ulong a1, a0, b1, b0;

        smul_ppmm(a1, a0, *p, *s);
        smul_ppmm(b1, b0, *q, *r);
        sub_ddmmss(a1, a0, a1, a0, b1, b0);

        if ((slong) a1 < 0)
            return -1;
        if ((slong) a1 > 0)
            return 1;
        return a0 != 0;
    }

    if (fmpz_equal(q, s))
        return fmpz_cmp(p, r);

    s1 = fmpz_sgn(p);
    s2 = fmpz_sgn(r);

    if (s1 != s2)
        return s1 < s2 ? -1 : 1;

    if (s1 == 0)
        return -s2;

    if (s2 == 0)
        return -s1;

    bp = fmpz_bits(p);
    bq = fmpz_bits(q);
    br = fmpz_bits(r);
    bs = fmpz_bits(s);

    if (bp + bs + 1 < br + bq)
        return -s1;

    if (bp + bs > br + bq + 1)
        return s1;

    if (fmpz_is_one(q))
    {
        fmpz_init(t);
        fmpz_mul(t, p, s);
        res = fmpz_cmp(t, r);
        fmpz_clear(t);
    }
    else if (fmpz_is_one(s))
    {
        fmpz_init(u);
        fmpz_mul(u, q, r);

        res = fmpz_cmp(p, u);

        fmpz_clear(u);
    }
    else
    {
        fmpz_init(t);
        fmpz_init(u);

        fmpz_mul(t, p, s);
        fmpz_mul(u, q, r);

        res = fmpz_cmp(t, u);

        fmpz_clear(t);
        fmpz_clear(u);
    }

    return res;
}

int
fmpq_cmp2(const fmpq_t x, const fmpq_t y)
{
    return _fmpq_cmp2(fmpq_numref(x), fmpq_denref(x),
                     fmpq_numref(y), fmpq_denref(y));
}

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("contains....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 100000 * arb_test_multiplier(); iter++)
    {
        arb_t a, b;
        fmpq_t am, ar, bm, br, t, u;
        int c1, c2;

        arb_init(a);
        arb_init(b);

        fmpq_init(am);
        fmpq_init(ar);
        fmpq_init(bm);
        fmpq_init(br);
        fmpq_init(t);
        fmpq_init(u);

        arb_randtest(a, state, 1 + n_randint(state, 500), 14);
        arb_randtest(b, state, 1 + n_randint(state, 500), 14);

        arf_get_fmpq(am, arb_midref(a));
        mag_get_fmpq(ar, arb_radref(a));
        arf_get_fmpq(bm, arb_midref(b));
        mag_get_fmpq(br, arb_radref(b));

        c1 = arb_contains(a, b);

        fmpq_sub(t, am, ar);
        fmpq_sub(u, bm, br);
        c2 = fmpq_cmp2(t, u) <= 0;

        fmpq_add(t, am, ar);
        fmpq_add(u, bm, br);
        c2 = c2 && (fmpq_cmp2(t, u) >= 0);

        if (c1 != c2)
        {
            flint_printf("FAIL:\n\n");
            flint_printf("a = "); arb_print(a); flint_printf("\n\n");
            flint_printf("b = "); arb_print(b); flint_printf("\n\n");
            flint_printf("am = "); fmpq_print(am); flint_printf("\n\n");
            flint_printf("ar = "); fmpq_print(ar); flint_printf("\n\n");
            flint_printf("bm = "); fmpq_print(bm); flint_printf("\n\n");
            flint_printf("br = "); fmpq_print(br); flint_printf("\n\n");
            flint_printf("t = "); fmpq_print(t); flint_printf("\n\n");
            flint_printf("u = "); fmpq_print(u); flint_printf("\n\n");
            flint_printf("c1 = %d, c2 = %d\n\n", c1, c2);
            flint_abort();
        }

        arb_clear(a);
        arb_clear(b);

        fmpq_clear(am);
        fmpq_clear(ar);
        fmpq_clear(bm);
        fmpq_clear(br);
        fmpq_clear(t);
        fmpq_clear(u);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
