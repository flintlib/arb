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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "gamma.h"

/* evaluate sum_{k=0}^{n-1} c_{k+1} x^k
  precision estimate assumes x <= 1/2 */
void
gamma_taylor_eval_fmprb(fmprb_t y, const fmprb_t x, long prec)
{
    long i, n, wp;
    fmprb_t t, u, v;
    fmpr_t z;

    n = gamma_taylor_coeffs_for_prec(prec);
    gamma_taylor_precompute(n, prec);

    fmprb_init(t);
    fmprb_init(u);
    fmprb_init(v);

    /* note: c_n is actually stored as coefficient n - 1 */
    fmprb_set(u, gamma_taylor_coeffs + n - 1);

    for (i = n - 1; i > 0; i--)
    {
        wp = prec + (gamma_taylor_bound_cached(i) - i);
        wp = FLINT_MIN(wp, prec + 10);
        wp = FLINT_MAX(wp, 10);

        if (prec > 2000)
        {
            fmprb_set_round(v, x, wp);
            fmprb_mul(u, u, v, wp);
        }
        else
        {
            fmprb_mul(u, u, x, wp);
        }

        fmprb_set_round(t, gamma_taylor_coeffs + i - 1, wp);
        fmprb_add(u, u, t, wp);
    }

    /* add error */
    fmpr_init(z);
    fmprb_get_abs_ubound_fmpr(z, x, FMPRB_RAD_PREC);
    gamma_taylor_bound_remainder(z, z, n);
    fmprb_add_error_fmpr(u, z);
    fmpr_clear(z);

    fmprb_set_round(y, u, prec);

    fmprb_clear(t);
    fmprb_clear(u);
    fmprb_clear(v);
}

