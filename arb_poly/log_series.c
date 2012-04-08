/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "arb_poly.h"

void
arb_poly_log_series(arb_poly_t z, const arb_poly_t x, long n)
{
    arb_poly_t t, u;

    arb_poly_init(t, arb_poly_prec(z));
    arb_poly_init(u, arb_poly_prec(z));

    arb_poly_inv_series(t, x, n);
    arb_poly_derivative(u, x);
    arb_poly_mullow(t, t, u, n - 1);
    arb_poly_integral(z, t);

    arb_poly_clear(t);
    arb_poly_clear(u);
}
