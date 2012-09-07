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

#include <math.h>
#include "arith.h"
#include "fmprb.h"

static void
fmprb_zeta_ui_mpfr(fmprb_t x, ulong n, long prec)
{
    mpfr_t t;
    mpfr_init2(t, prec);
    mpfr_zeta_ui(t, n, MPFR_RNDD);
    fmprb_zero(x);
    fmpr_set_mpfr(fmprb_midref(x), t);
    fmprb_add_error_2exp_si(x, mpfr_get_exp(t)-prec);
    mpfr_clear(t);
}

void
fmprb_zeta_ui(fmprb_t x, ulong n, long prec)
{
    /* asymptotic case */
    if (n == 0 || n > 0.7 * prec)
    {
        fmprb_zeta_ui_mpfr(x, n, prec);
    }
    else if (n == 3)
    {
        fmprb_const_zeta3_bsplit(x, prec);
    }
    /* small even n */
    else if ((n % 2 == 0) && (n < 40 + 0.11*prec))
    {
        fmprb_zeta_ui_bernoulli(x, n, prec);
    }
    /* small odd n, extremely high precision */
    else if (n < prec * 0.0006)
    {
        fmprb_zeta_ui_bsplit(x, n, prec);
    }
    /* large n */
    else if (prec > 20 && n > 0.4 * pow(prec, 0.8))
    {
        fmprb_zeta_ui_euler_product(x, n, prec);
    }
    /* fallback */
    else
    {
        fmprb_zeta_ui_mpfr(x, n, prec);
    }
}
