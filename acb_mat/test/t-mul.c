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

#include "acb_mat.h"


ACB_INLINE int
_acb_contains_fmpq_fmpq(const acb_t x, const fmpq_t a, const fmpq_t b)
{
    return arb_contains_fmpq(acb_realref(x), a) &&
            arb_contains_fmpq(acb_imagref(x), b);
}

ACB_INLINE void
_acb_set_fmpq_fmpq(acb_t z, const fmpq_t x, const fmpq_t y, slong prec)
{
    arb_set_fmpq(acb_realref(z), x, prec);
    arb_set_fmpq(acb_imagref(z), y, prec);
}

int
_acb_mat_contains_fmpq_mat_fmpq_mat(const acb_mat_t mat1,
    const fmpq_mat_t mat2, const fmpq_mat_t mat3)
{
    slong i, j, n, m;

    n = acb_mat_nrows(mat1);
    m = acb_mat_ncols(mat1);

    if ((n != acb_mat_nrows(mat2)) || (m != acb_mat_ncols(mat2)) ||
        (n != acb_mat_nrows(mat3)) || (m != acb_mat_ncols(mat3)))
        return 0;

    for (i = 0; i < n; i++)
        for (j = 0; j < m; j++)
            if (!_acb_contains_fmpq_fmpq(acb_mat_entry(mat1, i, j),
                fmpq_mat_entry(mat2, i, j), fmpq_mat_entry(mat3, i, j)))
                return 0;

    return 1;
}

void
_acb_mat_set_fmpq_mat_fmpq_mat(acb_mat_t dest,
    const fmpq_mat_t realsrc, const fmpq_mat_t imagsrc, slong prec)
{
    slong i, j;

    if (acb_mat_ncols(dest) != 0)
    {
        for (i = 0; i < acb_mat_nrows(dest); i++)
            for (j = 0; j < acb_mat_ncols(dest); j++)
                _acb_set_fmpq_fmpq(acb_mat_entry(dest, i, j),
                    fmpq_mat_entry(realsrc, i, j),
                    fmpq_mat_entry(imagsrc, i, j), prec);
    }
}


int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("mul....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        slong m, n, k;
        slong qbits1, qbits2;
        slong rbits1, rbits2, rbits3;
        fmpq_mat_t A, B, C, D, E, F, T;
        acb_mat_t a, b, c, d;

        qbits1 = 2 + n_randint(state, 200);
        qbits2 = 2 + n_randint(state, 200);
        rbits1 = 2 + n_randint(state, 200);
        rbits2 = 2 + n_randint(state, 200);
        rbits3 = 2 + n_randint(state, 200);

        m = n_randint(state, 10);
        n = n_randint(state, 10);
        k = n_randint(state, 10);

        fmpq_mat_init(A, m, n);
        fmpq_mat_init(B, m, n);
        fmpq_mat_init(C, n, k);
        fmpq_mat_init(D, n, k);
        fmpq_mat_init(E, m, k);
        fmpq_mat_init(F, m, k);
        fmpq_mat_init(T, m, k);

        acb_mat_init(a, m, n);
        acb_mat_init(b, n, k);
        acb_mat_init(c, m, k);
        acb_mat_init(d, m, k);

        fmpq_mat_randtest(A, state, qbits1);
        fmpq_mat_randtest(B, state, qbits1);
        fmpq_mat_randtest(C, state, qbits2);
        fmpq_mat_randtest(D, state, qbits2);

        /* E, F = AC - BD, AD + BC */
        fmpq_mat_mul(E, A, C);
        fmpq_mat_mul(T, B, D);
        fmpq_mat_sub(E, E, T);
        fmpq_mat_mul(F, A, D);
        fmpq_mat_mul(T, B, C);
        fmpq_mat_add(F, F, T);

        /* a, b = A + B*i, C + D*i */
        _acb_mat_set_fmpq_mat_fmpq_mat(a, A, B, rbits1);
        _acb_mat_set_fmpq_mat_fmpq_mat(b, C, D, rbits2);
        acb_mat_mul(c, a, b, rbits3);

        if (!_acb_mat_contains_fmpq_mat_fmpq_mat(c, E, F))
        {
            flint_printf("FAIL\n\n");
            flint_printf("m = %wd, n = %wd, k = %wd, bits3 = %wd\n", m, n, k, rbits3);

            flint_printf("A = "); fmpq_mat_print(A); flint_printf("\n\n");
            flint_printf("B = "); fmpq_mat_print(B); flint_printf("\n\n");
            flint_printf("C = "); fmpq_mat_print(C); flint_printf("\n\n");
            flint_printf("D = "); fmpq_mat_print(D); flint_printf("\n\n");
            flint_printf("E = "); fmpq_mat_print(E); flint_printf("\n\n");
            flint_printf("F = "); fmpq_mat_print(F); flint_printf("\n\n");

            flint_printf("a = "); acb_mat_printd(a, 15); flint_printf("\n\n");
            flint_printf("b = "); acb_mat_printd(b, 15); flint_printf("\n\n");
            flint_printf("c = "); acb_mat_printd(c, 15); flint_printf("\n\n");

            abort();
        }

        /* test aliasing with a */
        if (acb_mat_nrows(a) == acb_mat_nrows(c) &&
            acb_mat_ncols(a) == acb_mat_ncols(c))
        {
            acb_mat_set(d, a);
            acb_mat_mul(d, d, b, rbits3);
            if (!acb_mat_equal(d, c))
            {
                flint_printf("FAIL (aliasing 1)\n\n");
                abort();
            }
        }

        /* test aliasing with b */
        if (acb_mat_nrows(b) == acb_mat_nrows(c) &&
            acb_mat_ncols(b) == acb_mat_ncols(c))
        {
            acb_mat_set(d, b);
            acb_mat_mul(d, a, d, rbits3);
            if (!acb_mat_equal(d, c))
            {
                flint_printf("FAIL (aliasing 2)\n\n");
                abort();
            }
        }

        fmpq_mat_clear(A);
        fmpq_mat_clear(B);
        fmpq_mat_clear(C);
        fmpq_mat_clear(D);
        fmpq_mat_clear(E);
        fmpq_mat_clear(F);
        fmpq_mat_clear(T);

        acb_mat_clear(a);
        acb_mat_clear(b);
        acb_mat_clear(c);
        acb_mat_clear(d);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
