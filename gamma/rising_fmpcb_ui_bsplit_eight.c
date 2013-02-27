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

/* x(x+1)...(x+7) = (28 + 98x + 63x^2 + 14x^3 + x^4)^2 - 16 (7+2x)^2 */
static void
rfac_eight(fmpcb_t t, const fmpcb_t x, long prec)
{
    fmpcb_t u, v;

    fmpcb_init(u);
    fmpcb_init(v);

    /* t = x^2, v = x^3, u = x^4 */
    fmpcb_mul(t, x, x, prec);
    fmpcb_mul(v, x, t, prec);
    fmpcb_mul(u, t, t, prec);

    /* u = (28 + ...)^2 */
    fmpcb_addmul_ui(u, v, 14UL, prec);
    fmpcb_addmul_ui(u, t, 63UL, prec);
    fmpcb_addmul_ui(u, x, 98UL, prec);
    fmpcb_add_ui(u, u, 28UL, prec);
    fmpcb_mul(u, u, u, prec);

    /* 16 (7+2x)^2 = 784 + 448x + 64x^2 */
    fmpcb_sub_ui(u, u, 784UL, prec);
    fmpcb_submul_ui(u, x, 448UL, prec);
    fmpcb_mul_2exp_si(t, t, 6);
    fmpcb_sub(t, u, t, prec);

    fmpcb_clear(u);
    fmpcb_clear(v);
}

/* assumes y and x not aliased, the length is a positive multiple of 8 */
static void
bsplit(fmpcb_t y, const fmpcb_t x, ulong a, ulong b, long prec)
{
    fmpcb_t t;
    fmpcb_init(t);

    if (b - a == 8)
    {
        fmpcb_add_ui(t, x, a, prec);
        rfac_eight(y, t, prec);
    }
    else
    {
        ulong m = a + ((b - a) / 16) * 8;
        bsplit(y, x, a, m, prec);
        bsplit(t, x, m, b, prec);
        fmpcb_mul(y, y, t, prec);
    }

    fmpcb_clear(t);
}

void
gamma_rising_fmpcb_ui_bsplit_eight(fmpcb_t y, const fmpcb_t x, ulong n, long prec)
{
    if (n == 0)
    {
        fmpcb_one(y);
    }
    else if (n == 1)
    {
        fmpcb_set_round(y, x, prec);
    }
    else
    {
        ulong k, a;
        long wp;
        fmpcb_t t, u;

        wp = FMPR_PREC_ADD(prec, FLINT_BIT_COUNT(n));

        fmpcb_init(t);
        fmpcb_init(u);

        if (n >= 8)
        {
            bsplit(t, x, 0, (n / 8) * 8, wp);
            a = (n / 8) * 8;
        }
        else
        {
            fmpcb_set(t, x);
            a = 1;
        }

        for (k = a; k < n; k++)
        {
            fmpcb_add_ui(u, x, k, wp);
            fmpcb_mul(t, t, u, wp);
        }

        fmpcb_set_round(y, t, prec);

        fmpcb_clear(t);
        fmpcb_clear(u);
    }
}

