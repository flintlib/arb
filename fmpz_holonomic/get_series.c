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

/* TODO: use the fmpz_poly functions */

void
_fmpz_poly_from_newton_basis(fmpz * poly, const fmpz * roots, long n)
{
    long i, j;

    for (i = n - 1; i > 0; i--)
        for (j = i - 1; j < n - 1; j++)
            fmpz_submul(poly + j, poly + j + 1, roots + i - 1);
}

void
_fmpz_poly_to_newton_basis(fmpz * poly, const fmpz * roots, long n)
{
    long i, j;

    for (i = 1; i < n; i++)
        for (j = n - 2; j > i - 2; j--)
            fmpz_addmul(poly + j, poly + j + 1, roots + i - 1);
}

void
fmpz_poly_from_newton_basis(fmpz_poly_t g, const fmpz_poly_t f, const
fmpz * roots)
{
    if (f != g)
        fmpz_poly_set(g, f);

    _fmpz_poly_from_newton_basis(g->coeffs, roots, g->length);
}

void
fmpz_poly_to_newton_basis(fmpz_poly_t g, const fmpz_poly_t f, const
fmpz * roots)
{
    if (f != g)
        fmpz_poly_set(g, f);

    _fmpz_poly_to_newton_basis(g->coeffs, roots, g->length);
}




void
fmpz_holonomic_get_series(fmpz_holonomic_t re, const fmpz_holonomic_t de)
{
    long i, j, k, start, d, len, r;
    fmpz * roots;

    /*
    Example:
    Input polynomials: a + b*D_x + c*D_x^2, degrees bounded by 4
    Output polynomials: r0, r1, r2, r3, r4, r5, r6

            Coefficients  .  Basis
     r0 =   a4               (1)
     r1 =   a3, b4           (1), (n+1)
     r2 =   a2, b3, c4       (1), (n+2), (n+2)(n+1)
     r3 =   a1, b2, c3       (1), (n+3), (n+3)(n+2)
     r4 =   a0, b1, c2       (1), (n+4), (n+4)(n+3)
     r5 =       b0, c1            (n+5), (n+5)(n+4)
     r6 =           c0                   (n+5)(n+6)

    Note that we can tell in advance if some low-order polynomials in the
    output will be zero. We denote this offset by start and include it
    in the rising factorial basis, to avoid the need for a final adjustment.
    */

    r = fmpz_holonomic_order(de);
    d = fmpz_holonomic_degree(de);
    len = r + 1;

    start = d + 1;
    for (k = 0; k < len; k++)
        start = FLINT_MIN(start, d - fmpz_poly_degree(de->coeffs + k) + k);

    roots = malloc(sizeof(fmpz) * (len - 1));

    fmpz_holonomic_fit_length(re, d + len - start);

    for (k = start; k < d + len; k++)
    {
        i = k - start;
        fmpz_poly_zero(re->coeffs + i);

        for (j = 0; j < len; j++)
        {
            long v = d + j - k;

            if (v >= 0 && v < fmpz_poly_length(de->coeffs + j))
                fmpz_poly_set_coeff_fmpz(re->coeffs + i, j,
                    (de->coeffs + j)->coeffs + v);
        }

        if (!fmpz_poly_is_zero(re->coeffs + i))
        {
            for (j = 0; j < len - 1; j++)
                roots[j] = -(i - j);
            fmpz_poly_from_newton_basis(re->coeffs + i, re->coeffs + i, roots);
        }
    }

    _fmpz_holonomic_set_length(re, d + len - start);
    fmpz_holonomic_seq_normalise(re);
    free(roots);
}

