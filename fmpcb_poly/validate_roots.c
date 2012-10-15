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

#include "fmpcb_poly.h"

/*
Given a list of approximate roots of the input polynomial, this
function sets a rigorous bounding interval for each root, and determines
which roots are isolated from all the other roots.
It then rearranges the list of roots so that the isolated roots
are at the front of the list, and returns the count of isolated roots.

If the return value equals the degree of the polynomial, then all
roots have been separated. If the return value is smaller, all the
remaining intervals are guaranteed to contain roots, but
it is possible that not all of the polynomial's roots are contained
among them.
*/

long
_fmpcb_poly_validate_roots(fmpcb_struct * roots,
        const fmpcb_struct * poly, long len, long prec)
{
    long i, j, deg;
    long isolated, nonisolated, total_isolated;
    fmpcb_struct * deriv;
    fmpcb_struct * tmp;
    int *overlap;

    deg = len - 1;

    deriv = _fmpcb_vec_init(deg);
    overlap = flint_calloc(deg, sizeof(int));
    tmp = flint_malloc(sizeof(fmpcb_struct) * deg);

    _fmpcb_poly_derivative(deriv, poly, len, prec);

    /* compute an inclusion interval for each point */
    for (i = 0; i < deg; i++)
    {
        _fmpcb_poly_root_inclusion(roots + i, roots + i,
            poly, deriv, len, prec);
    }

    /* find which points do not overlap with any other points */
    for (i = 0; i < deg; i++)
    {
        for (j = i + 1; j < deg; j++)
        {
            if (fmpcb_overlaps(roots + i, roots + j))
            {
                overlap[i] = overlap[j] = 1;
            }
        }
    }

    /* count and move all isolated roots to the front of the array */
    total_isolated = 0;
    for (i = 0; i < deg; i++)
        total_isolated += (overlap[i] == 0);

    for (i = 0; i < deg; i++)
        tmp[i] = roots[i];

    isolated = 0;
    nonisolated = 0;
    for (i = 0; i < deg; i++)
    {
        if (overlap[i] == 0)
        {
            roots[isolated] = tmp[i];
            isolated++;
        }
        else
        {
            roots[total_isolated + nonisolated] = tmp[i];
            nonisolated++;
        }
    }

    _fmpcb_vec_clear(deriv, deg);
    flint_free(tmp);
    flint_free(overlap);

    return isolated;
}

