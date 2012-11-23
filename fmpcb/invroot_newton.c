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

#include "fmpcb.h"

/*
Given one inverse mth root r0 (with a valid error bound) of the complex
number a, assumed to be isolated from the conjugate roots and from the origin,
lift it from precision startprec to prec using Newton iteration,
solving f(z) = (1/z)^m - a = 0.

Given an error bound e_n for an input term z_n at step n,
we have the propagated output error
e_{n+1} < |1/f'(z)| * sum_{k>=2} |f^{(k)}| / k! * (e_n)^k.
Replacing k! by (k-2)! gives the closed form used below.

*/

void
fmpcb_invroot_newton(fmpcb_t r, const fmpcb_t a, ulong m,
    const fmpcb_t r0, long startprec, long prec)
{
    long iters, precs[FLINT_BITS];
    long i, extra, wp, rad_prec;

    fmpr_t en, enew, zlo, zhi, v;
    fmpcb_t t, z;

    fmpr_init(en);
    fmpr_init(enew);
    fmpr_init(zlo);
    fmpr_init(zhi);
    fmpr_init(v);

    fmpcb_init(t);
    fmpcb_init(z);

    fmpcb_set(z, r0);

    rad_prec = FMPRB_RAD_PREC;

    extra = 10 + FLINT_BIT_COUNT(m);
    precs[0] = prec + extra;
    iters = 1;

    while ((iters < FLINT_BITS) && (precs[iters-1] + extra > 2*startprec))
    {
        precs[iters] = (precs[iters-1] / 2) + extra;
        iters++;
    }

    for (i = iters - 1; i >= 0; i--)
    {
        wp = precs[i];

        /* printf("lifting to precision %ld\n", wp); */

        /* bounds for old error */
        fmpcb_get_rad_ubound_fmpr(en, z, rad_prec);

        fmpcb_get_abs_lbound_fmpr(zlo, z, rad_prec);
        fmpcb_get_abs_ubound_fmpr(zhi, z, rad_prec);

        /* to improve, we require |z| > 0, i.e. en < |z| */
        if (fmpr_cmp(en, zlo) >= 0)
            break;

        /* z[n+1] = z[n] * (m + 1 - a * z[n]^m) / m */
        fmpr_zero(fmprb_radref(fmpcb_realref(z)));
        fmpr_zero(fmprb_radref(fmpcb_imagref(z)));

        fmpcb_pow_ui(t, z, m, wp);
        fmpcb_mul(t, t, a, wp);
        fmpcb_neg(t, t);
        fmpcb_add_ui(t, t, m + 1, wp);
        fmpcb_mul(t, t, z, wp);
        fmpcb_div_ui(t, t, m, wp);

        /* new error bound = en^2 * (m+1) / (|zn| / (|zn| - en))^(m+2) / |zn| */
        fmpr_mul(enew, en, en, rad_prec, FMPR_RND_UP);
        fmpr_mul_ui(enew, enew, m + 1, rad_prec, FMPR_RND_UP);
        fmpr_mul(enew, enew, zhi, rad_prec, FMPR_RND_UP);
        fmpr_sub(v, zlo, en, rad_prec, FMPR_RND_DOWN);
        fmpr_div(v, zhi, v, rad_prec, FMPR_RND_UP);
        fmpr_pow_sloppy_ui(v, v, m + 2, rad_prec, FMPR_RND_UP);
        fmpr_div(v, v, zlo, rad_prec, FMPR_RND_UP);
        fmpr_mul(enew, enew, v, rad_prec, FMPR_RND_UP);

        /* quit if there was no improvement */
        if (fmpr_cmp(enew, en) >= 0)
            break;

        fmpcb_set(z, t);

        /* add upper bound for new radius */
        fmprb_add_error_fmpr(fmpcb_realref(z), enew);
        fmprb_add_error_fmpr(fmpcb_imagref(z), enew);
    }

    fmpcb_set(r, z);

    fmpcb_clear(t);
    fmpcb_clear(z);

    fmpr_clear(en);
    fmpr_clear(enew);
    fmpr_clear(zlo);
    fmpr_clear(zhi);
    fmpr_clear(v);
}

