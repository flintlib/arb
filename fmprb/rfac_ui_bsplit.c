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

#include "fmprb.h"

/* x(x+1)...(x+7) = (28 + 98x + 63x^2 + 14x^3 + x^4)^2 - 16 (7+2x)^2 */
static void
rfac_eight(fmprb_t t, const fmprb_t x, long prec)
{
    fmprb_t u, v;

    fmprb_init(u);
    fmprb_init(v);

    fmprb_mul(t, x, x, prec);
    fmprb_mul(v, x, t, prec);
    fmprb_mul(u, t, t, prec);

    fmprb_addmul_ui(u, v, 14UL, prec);
    fmprb_addmul_ui(u, t, 63UL, prec);
    fmprb_addmul_ui(u, x, 98UL, prec);
    fmprb_add_ui(u, u, 28UL, prec);

    fmprb_mul(u, u, u, prec);

    fmprb_mul_2exp_si(t, x, 1);
    fmprb_add_ui(t, t, 7UL, prec);
    fmprb_mul(t, t, t, prec);
    fmprb_mul_2exp_si(t, t, 4);

    fmprb_sub(t, u, t, prec);

    fmprb_clear(u);
    fmprb_clear(v);
}

/* assumes that the length is a multiple of 8 */
static void
bsplit_eight(fmprb_t y, const fmprb_t x,
    ulong a, ulong b, long prec)
{
    if (b - a == 8)
    {
        fmprb_t t;
        fmprb_init(t);
        fmprb_add_ui(t, x, a, prec);
        rfac_eight(y, t, prec);
        fmprb_clear(t);
    }
    else
    {
        ulong m = a + ((b - a) / 16) * 8;

        fmprb_t L, R;

        fmprb_init(L);
        fmprb_init(R);

        bsplit_eight(L, x, a, m, prec);
        bsplit_eight(R, x, m, b, prec);
        fmprb_mul(y, L, R, prec);

        fmprb_clear(L);
        fmprb_clear(R);
    }
}

void
fmprb_rfac_ui_bsplit(fmprb_t y, const fmprb_t x, ulong n, long prec)
{
    if (n == 0)
    {
        fmprb_one(y);
    }
    else if (n == 1)
    {
        fmprb_set(y, x);
    }
    else
    {
        long k, a, wp;
        fmprb_t t, u;

        wp = FMPR_PREC_ADD(prec, FLINT_BIT_COUNT(n));

        fmprb_init(t);
        fmprb_init(u);

        if (n >= 8)
        {
            bsplit_eight(t, x, 0, (n / 8) * 8, wp);
            a = (n / 8) * 8;
        }
        else
        {
            fmprb_set(t, x);
            a = 1;
        }

        for (k = a; k < n; k++)
        {
            fmprb_add_ui(u, x, k, wp);
            fmprb_mul(t, t, u, wp);
        }

        fmprb_set(y, t);

        fmprb_clear(t);
        fmprb_clear(u);
    }
}
