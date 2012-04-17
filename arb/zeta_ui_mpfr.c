/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "arb.h"

void
arb_zeta_ui_mpfr(arb_t x, ulong n)
{
    mpfr_t t;

    mpfr_init2(t, FLINT_MAX(arb_prec(x), FLINT_BITS));

    if (n == 1)
        mpfr_const_euler(t, MPFR_RNDN);
    else
        mpfr_zeta_ui(t, n, MPFR_RNDN);

    arb_set_mpfr(x, t, 1);
    mpfr_clear(t);
}
