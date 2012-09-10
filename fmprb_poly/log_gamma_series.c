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

#include "fmprb_poly.h"


/* series expansion for log(gamma(1-x)) at x = 0 */

void
fmprb_poly_log_gamma_series(fmprb_poly_t z, long n, long prec)
{
    long i;

    fmprb_poly_fit_length(z, n);
    _fmprb_poly_set_length(z, n);

    if (n > 0) fmprb_zero(z->coeffs);
    if (n > 1) fmprb_const_euler_brent_mcmillan(z->coeffs + 1, prec);
    if (n > 2) fmprb_zeta_ui_vec(z->coeffs + 2, 2, n - 2, prec);

    for (i = 2; i < n; i++)
        fmprb_div_ui(z->coeffs + i, z->coeffs + i, i, prec);

    _fmprb_poly_normalise(z);
}
