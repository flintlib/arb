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

#include "acb_poly.h"

long
_acb_poly_validate_roots(acb_ptr roots,
        acb_srcptr poly, long len, long prec)
{
    long i, j, deg;
    long isolated, nonisolated, total_isolated;
    acb_ptr deriv;
    acb_ptr tmp;
    int *overlap;

    deg = len - 1;

    deriv = _acb_vec_init(deg);
    overlap = flint_calloc(deg, sizeof(int));
    tmp = flint_malloc(sizeof(acb_struct) * deg);

    _acb_poly_derivative(deriv, poly, len, prec);

    /* compute an inclusion interval for each point */
    for (i = 0; i < deg; i++)
    {
        _acb_poly_root_inclusion(roots + i, roots + i,
            poly, deriv, len, prec);
    }

    /* find which points do not overlap with any other points */
    for (i = 0; i < deg; i++)
    {
        for (j = i + 1; j < deg; j++)
        {
            if (acb_overlaps(roots + i, roots + j))
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

    _acb_vec_clear(deriv, deg);
    flint_free(tmp);
    flint_free(overlap);

    return isolated;
}

