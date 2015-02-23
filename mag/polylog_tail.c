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

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include "mag.h"

void
mag_polylog_tail(mag_t u, const mag_t z, long sigma, ulong d, ulong N)
{
    mag_t TN, UN, t;

    if (N < 2)
    {
        mag_inf(u);
        return;
    }

    mag_init(TN);
    mag_init(UN);
    mag_init(t);

    if (mag_cmp_2exp_si(z, 0) >= 0)
    {
        mag_inf(u);
    }
    else
    {
        /* Bound T(N) */
        mag_pow_ui(TN, z, N);

        /* multiply by log(N)^d */
        if (d > 0)
        {
            mag_log_ui(t, N);
            mag_pow_ui(t, t, d);
            mag_mul(TN, TN, t);
        }

        /* multiply by 1/k^s */
        if (sigma > 0)
        {
            mag_set_ui_lower(t, N);
            mag_pow_ui_lower(t, t, sigma);
            mag_div(TN, TN, t);
        }
        else if (sigma < 0)
        {
            mag_set_ui(t, N);
            mag_pow_ui(t, t, -sigma);
            mag_mul(TN, TN, t);
        }

        /* Bound U(N) */
        mag_set(UN, z);

        /* multiply by (1 + 1/N)**S */
        if (sigma < 0)
        {
            mag_binpow_uiui(t, N, -sigma);
            mag_mul(UN, UN, t);
        }

        /* multiply by (1 + 1/(N log(N)))^d */
        if (d > 0)
        {
            ulong nl;

            /* rounds down */
            nl = mag_d_log_lower_bound(N) * N * (1 - 1e-13);

            mag_binpow_uiui(t, nl, d);
            mag_mul(UN, UN, t);
        }

        /* T(N) / (1 - U(N)) */
        if (mag_cmp_2exp_si(UN, 0) >= 0)
        {
            mag_inf(u);
        }
        else
        {
            mag_one(t);
            mag_sub_lower(t, t, UN);
            mag_div(u, TN, t);
        }
    }

    mag_clear(TN);
    mag_clear(UN);
    mag_clear(t);
}

