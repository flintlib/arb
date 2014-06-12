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

#include "hypgeom.h"

void
hypgeom_init(hypgeom_t hyp)
{
    fmpz_poly_init(hyp->A);
    fmpz_poly_init(hyp->B);
    fmpz_poly_init(hyp->P);
    fmpz_poly_init(hyp->Q);
    mag_init(hyp->MK);
    hyp->have_precomputed = 0;
}

void
hypgeom_clear(hypgeom_t hyp)
{
    fmpz_poly_clear(hyp->A);
    fmpz_poly_clear(hyp->B);
    fmpz_poly_clear(hyp->P);
    fmpz_poly_clear(hyp->Q);
    mag_clear(hyp->MK);
}
