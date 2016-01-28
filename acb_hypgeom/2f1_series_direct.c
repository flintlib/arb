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

#include "acb_hypgeom.h"

void
acb_hypgeom_2f1_series_direct(acb_poly_t res,
    const acb_poly_t a, const acb_poly_t b,
    const acb_poly_t c, const acb_poly_t z, int regularized, slong len, slong prec)
{
    acb_poly_struct aa[4];

    aa[0] = *a;
    aa[1] = *b;
    aa[2] = *c;
    acb_poly_init(&aa[3]);
    acb_poly_one(&aa[3]);

    acb_hypgeom_pfq_series_direct(res, aa, 2, aa + 2, 2, z,
        regularized, -1, len, prec);

    acb_poly_clear(&aa[3]);
}

