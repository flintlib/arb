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
#include "fmprb_poly.h"
#include "fmpcb_poly.h"

/*
    4 / (2pi)^(2M) * (N+a)^(1-s-2M) * rf(s, 2+2M-1) * (N+a)^(-x)
    Note: this is with the tail sum done up to M inclusive.

    For convergence, Re(N+a) > 0 and Re(1-s-2M) < 0.
 */


void
fmpcb_zeta_series_em_vec_bound(fmprb_struct * vec,
        const fmpcb_t s, const fmpcb_t a, long d, long N, long M, long wp)
{
    fmprb_struct *t, *u;
    fmprb_struct sx[2];
    fmprb_t x, y, z;
    fmpcb_t w, h, Na, Ms;
    long i;

    fmpcb_init(Na);
    fmpcb_init(Ms);

    fmpcb_add_ui(Na, a, N, wp);
    fmpcb_add_ui(Ms, s, 2*M-1, wp);
    fmpcb_neg(Ms, Ms);

    if (!fmprb_is_positive(fmpcb_realref(Na)) ||
        !fmprb_is_negative(fmpcb_realref(Ms)) || N < 1 || M < 1)
    {
        fmpcb_clear(Na);
        fmpcb_clear(Ms);

        for (i = 0; i < d; i++)
        {
            fmpr_pos_inf(fmprb_midref(vec + i));
            fmpr_zero(fmprb_radref(vec + i));
        }
        return;
    }

    t = _fmprb_vec_init(d);
    u = _fmprb_vec_init(d);

    fmprb_init(x);
    fmprb_init(y);
    fmprb_init(z);

    fmpcb_init(w);
    fmpcb_init(h);

    /* x = 4 / (2pi)^(2M) */
    fmprb_const_pi(x, wp);
    fmprb_mul_ui(x, x, 2, wp);
    fmprb_pow_ui(x, x, 2 * M, wp);
    fmprb_ui_div(x, 4, x, wp);

    /* y = |(N+a)^(-(s+2M-1))| */
    fmpcb_pow(w, Na, Ms, wp);
    fmpcb_abs(y, w, wp);

    /* x = x * y */
    fmprb_mul(x, x, y, wp);

    /* rising factorial */
    fmprb_init(sx + 0);
    fmprb_init(sx + 1);
    fmpcb_abs(sx + 0, s, wp);
    fmprb_one(sx + 1);
    _fmprb_poly_rfac_series_ui(t, sx, 2, 2+2*M-1, d, wp);
    fmprb_clear(sx + 0);
    fmprb_clear(sx + 1);

    /* exponential */
    fmpcb_log(w, Na, wp);
    fmpcb_abs(z, w, wp);
    fmprb_one(u + 0);
    for (i = 1; i < d; i++)
    {
        fmprb_mul(u + i, u + i - 1, z, wp);
        fmprb_div_ui(u + i, u + i, i, wp);
    }

    /* multiply everything together */
    _fmprb_poly_mullow(vec, t, d, u, d, d, wp);
    _fmprb_vec_scalar_mul(vec, vec, d, x, wp);

    fmprb_clear(x);
    fmprb_clear(y);
    fmprb_clear(z);

    fmpcb_clear(w);
    fmpcb_clear(h);

    fmpcb_clear(Na);
    fmpcb_clear(Ms);

    _fmprb_vec_clear(t, d);
    _fmprb_vec_clear(u, d);
}

void
fmpcb_zeta_series_em_bound(fmpr_t bound,
        const fmpcb_t s, const fmpcb_t a, long d, long N, long M, long wp)
{
    fmprb_struct * vec = _fmprb_vec_init(d);
    fmpcb_zeta_series_em_vec_bound(vec, s, a, d, N, M, wp);
    _fmprb_vec_get_abs_ubound_fmpr(bound, vec, d, wp);
    _fmprb_vec_clear(vec, d);
}

