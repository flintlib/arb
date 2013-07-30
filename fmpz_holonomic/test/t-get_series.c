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

#include "fmpz_holonomic.h"

void
check(const fmpz_holonomic_t re,
    const fmpz_holonomic_t ans, const fmpz_holonomic_t de)
{
    if (!fmpz_holonomic_equal(re, ans))
    {
        printf("FAIL: ");
        printf("de: "); fmpz_holonomic_print(de, "x", "Dx"); printf("\n\n");
        printf("re: "); fmpz_holonomic_print(re, "n", "Sn"); printf("\n\n");
        printf("ans: "); fmpz_holonomic_print(ans, "n", "Sn"); printf("\n\n");
        abort();
    }
}

int main()
{
    printf("get_series....");
    fflush(stdout);

    {
        fmpq_t s;
        fmpz_holonomic_t re, de, ans;

        fmpq_init(s);
        fmpz_holonomic_init(re);
        fmpz_holonomic_init(de);
        fmpz_holonomic_init(ans);
        fmpz_holonomic_fit_length(ans, 4);

        fmpq_set_si(s, -7, 3);

        /* exp(x) */
        fmpz_holonomic_fun_set_exp(de);
        fmpz_holonomic_get_series(re, de);
        fmpz_poly_set_si(ans->coeffs + 0, -1);
        fmpz_poly_set_si2(ans->coeffs + 1, 1, 1);
        _fmpz_holonomic_set_length(ans, 2);
        check(re, ans, de);

        /* sin/cos(x) */
        fmpz_holonomic_fun_set_sin_cos(de);
        fmpz_holonomic_get_series(re, de);
        fmpz_poly_set_str(ans->coeffs + 0, "1  1");
        fmpz_poly_set_str(ans->coeffs + 1, "0");
        fmpz_poly_set_str(ans->coeffs + 2, "3  2 3 1");
        _fmpz_holonomic_set_length(ans, 3);
        check(re, ans, de);

        /* log(x) */
        fmpz_holonomic_fun_set_log(de);
        fmpz_holonomic_get_series(re, de);
        fmpz_poly_set_si(ans->coeffs + 0, 1);
        _fmpz_holonomic_set_length(ans, 1);
        check(re, ans, de);

        /* log(s+x) */
        fmpz_holonomic_shift_fmpq(de, de, s);
        fmpz_holonomic_get_series(re, de);
        fmpz_poly_set_str(ans->coeffs + 0, "2  0 -3");
        fmpz_poly_set_str(ans->coeffs + 1, "2  7 7");
        _fmpz_holonomic_set_length(ans, 2);
        check(re, ans, de);

        /* atan(x) */
        fmpz_holonomic_fun_set_atan(de);
        fmpz_holonomic_get_series(re, de);
        fmpz_poly_set_str(ans->coeffs + 0, "2  0 1");
        fmpz_poly_set_str(ans->coeffs + 1, "0");
        fmpz_poly_set_str(ans->coeffs + 2, "2  2 1");
        _fmpz_holonomic_set_length(ans, 3);
        check(re, ans, de);

        /* atan(s+x) */
        fmpz_holonomic_shift_fmpq(de, de, s);
        fmpz_holonomic_get_series(re, de);
        fmpz_poly_set_str(ans->coeffs + 0, "2  0 9");
        fmpz_poly_set_str(ans->coeffs + 1, "2  -42 -42");
        fmpz_poly_set_str(ans->coeffs + 2, "2  116 58");
        _fmpz_holonomic_set_length(ans, 3);
        check(re, ans, de);

        /* erf(x) */
        fmpz_holonomic_fun_set_erf(de);
        fmpz_holonomic_get_series(re, de);
        fmpz_poly_set_str(ans->coeffs + 0, "2  0 2");
        fmpz_poly_set_str(ans->coeffs + 1, "0");
        fmpz_poly_set_str(ans->coeffs + 2, "3  2 3 1");
        _fmpz_holonomic_set_length(ans, 3);
        check(re, ans, de);

        /* erf(s+x) */
        fmpz_holonomic_shift_fmpq(de, de, s);
        fmpz_holonomic_get_series(re, de);
        fmpz_poly_set_str(ans->coeffs + 0, "2  0 6");
        fmpz_poly_set_str(ans->coeffs + 1, "2  -14 -14");
        fmpz_poly_set_str(ans->coeffs + 2, "3  6 9 3");
        _fmpz_holonomic_set_length(ans, 3);
        check(re, ans, de);

        /* x^s */
        fmpz_holonomic_fun_set_pow_fmpq(de, s);
        fmpz_holonomic_get_series(re, de);
        fmpz_poly_set_si(ans->coeffs + 0, 1);
        _fmpz_holonomic_set_length(ans, 1);
        check(re, ans, de);

        /* (x+s)^s */
        fmpz_holonomic_shift_fmpq(de, de, s);
        fmpz_holonomic_normalise_content(de);
        fmpz_holonomic_get_series(re, de);
        fmpz_poly_set_str(ans->coeffs + 0, "2  -7 -3");
        fmpz_poly_set_str(ans->coeffs + 1, "2  7 7");
        _fmpz_holonomic_set_length(ans, 2);
        check(re, ans, de);

        /* x^n */
        fmpz_holonomic_fun_set_pow_fmpz(de, fmpq_numref(s));
        fmpz_holonomic_get_series(re, de);
        fmpz_poly_set_si(ans->coeffs + 0, 1);
        _fmpz_holonomic_set_length(ans, 1);
        check(re, ans, de);

        /* (x+s)^n */
        fmpz_holonomic_shift_fmpq(de, de, s);
        fmpz_holonomic_normalise_content(de);
        fmpz_holonomic_get_series(re, de);
        fmpz_poly_set_str(ans->coeffs + 0, "2  -21 -3");
        fmpz_poly_set_str(ans->coeffs + 1, "2  7 7");
        _fmpz_holonomic_set_length(ans, 2);
        check(re, ans, de);

        fmpz_holonomic_clear(re);
        fmpz_holonomic_clear(de);
        fmpz_holonomic_clear(ans);
        fmpq_clear(s);
    }

    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

