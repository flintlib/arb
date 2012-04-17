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

/* series expansion for log(gamma(1-x)) at x = 0 */

void
arb_poly_gamma_log_series(arb_poly_t z, ulong n)
{
    arb_t c;
    long i, prec;

    prec = arb_poly_prec(z) + FLINT_BIT_COUNT(n);

    arb_init(c, prec);

    arb_poly_zero(z);
    _arb_poly_fit_length(z, n);
    z->length = n;

    fmpz_set_si(arb_poly_expref(z), -prec);

    for (i = 0; i < n; i++)
    {
        if (i == 0)
        {
            arb_zero(c);
        }
        else if (i == 1)
        {
            arb_const_euler_brent_mcmillan(c);
        }
        else
        {
            /* printf("%ld of %ld\n", i, n); */
            arb_zeta_ui(c, i);
            fmpz_tdiv_q_ui(arb_midref(c), arb_midref(c), i);
            fmpz_add_ui(arb_radref(c), arb_radref(c), 1);
        }

        _arb_poly_set_coeff_same_exp(z, i, c);
    }

    arb_clear(c);
}
