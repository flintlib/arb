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

static __inline__ void
fmprb_div_zero(fmprb_t z)
{
    fmpr_zero(fmprb_midref(z));
    fmpr_pos_inf(fmprb_radref(z));
    return;
}

void
fmprb_div_fmpr(fmprb_t z, const fmprb_t x, const fmpr_t y, long prec)
{
    long r;

    if (fmpr_is_zero(y))
    {
        fmprb_div_zero(z);
    }
    else if (fmprb_is_exact(x))
    {
        r = fmpr_div(fmprb_midref(z), fmprb_midref(x), y, prec, FMPR_RND_DOWN);
        fmpr_set_error_result(fmprb_radref(z), fmprb_midref(z), r);
    }
    else if (fmpr_is_pos_inf(fmprb_radref(x)))
    {
        fmpr_div(fmprb_midref(z), fmprb_midref(x), y, prec, FMPR_RND_DOWN);
        fmpr_pos_inf(fmprb_radref(z));
    }
    else
    {
        /* (x + a) / y = x/y + a/y */

        fmpr_divappr_abs_ubound(fmprb_radref(z), fmprb_radref(x), y, FMPRB_RAD_PREC);
        fmpr_abs(fmprb_radref(z), fmprb_radref(z));

        r = fmpr_div(fmprb_midref(z), fmprb_midref(x), y, prec, FMPR_RND_DOWN); 
        fmpr_add_error_result(fmprb_radref(z), fmprb_radref(z),
            fmprb_midref(z), r, FMPRB_RAD_PREC, FMPR_RND_UP);
    }
}

void
fmprb_div(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)
{
    long r;

    if (fmprb_is_exact(y))
    {
        fmprb_div_fmpr(z, x, fmprb_midref(y), prec);
    }
    else if (fmpr_is_pos_inf(fmprb_radref(x)) || fmpr_is_pos_inf(fmprb_radref(y)))
    {
        fmpr_div(fmprb_midref(z), fmprb_midref(x), fmprb_midref(y), prec, FMPR_RND_DOWN);
        fmpr_pos_inf(fmprb_radref(z));
    }
    else
    {
        fmpr_t t, u;

        fmpr_init(t);
        fmpr_init(u);

        /* numerator of error bound: |xb| + |ya|, rounded up */
        fmpr_mul(t, fmprb_midref(x), fmprb_radref(y), FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_abs(t, t);
        fmpr_mul(u, fmprb_midref(y), fmprb_radref(x), FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_abs(u, u);
        fmpr_add(t, t, u, FMPRB_RAD_PREC, FMPR_RND_UP);

        /* denominator of error bound: |y|(|y|-b), rounded down */
        if (fmpr_sgn(fmprb_midref(y)) > 0)
        {
            fmpr_sub(u, fmprb_midref(y), fmprb_radref(y), FMPRB_RAD_PREC, FMPR_RND_DOWN);
        }
        else
        {
            fmpr_add(u, fmprb_midref(y), fmprb_radref(y), FMPRB_RAD_PREC, FMPR_RND_DOWN);
            fmpr_neg(u, u);
        }

        if (fmpr_sgn(u) <= 0 || fmpr_is_nan(u))
        {
            fmprb_div_zero(z);
        }
        else
        {
            fmpr_mul(u, u, fmprb_midref(y), FMPRB_RAD_PREC, FMPR_RND_DOWN);
            fmpr_divappr_abs_ubound(t, t, u, FMPRB_RAD_PREC);

            r = fmpr_div(fmprb_midref(z), fmprb_midref(x), fmprb_midref(y), prec, FMPR_RND_DOWN);
            fmpr_add_error_result(fmprb_radref(z), t,
                fmprb_midref(z), r, FMPRB_RAD_PREC, FMPR_RND_UP);
        }

        fmpr_clear(t);
        fmpr_clear(u);
    }

    fmprb_adjust(z);
}

void
fmprb_div_ui(fmprb_t z, const fmprb_t x, ulong y, long prec)
{
    fmpr_t t;
    fmpr_init(t);
    fmpr_set_ui(t, y);
    fmprb_div_fmpr(z, x, t, prec);
    fmpr_clear(t);
}

void
fmprb_div_si(fmprb_t z, const fmprb_t x, long y, long prec)
{
    fmpr_t t;
    fmpr_init(t);
    fmpr_set_si(t, y);
    fmprb_div_fmpr(z, x, t, prec);
    fmpr_clear(t);
}

void
fmprb_div_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, long prec)
{
    fmpr_t t;
    fmpr_init(t);
    fmpr_set_fmpz(t, y);
    fmprb_div_fmpr(z, x, t, prec);
    fmpr_clear(t);
}

void
fmprb_fmpz_div_fmpz(fmprb_t y, const fmpz_t num, const fmpz_t den, long prec)
{
    fmpr_t p, q;
    long r;

    fmpr_init(p);
    fmpr_init(q);

    fmpr_set_fmpz(p, num);
    fmpr_set_fmpz(q, den);

    r = fmpr_div(fmprb_midref(y), p, q, prec, FMPR_RND_DOWN);
    fmpr_set_error_result(fmprb_radref(y), fmprb_midref(y), r);

    fmpr_clear(p);
    fmpr_clear(q);
}

void
fmprb_ui_div(fmprb_t z, ulong x, const fmprb_t y, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_set_ui(t, x);
    fmprb_div(z, t, y, prec);
    fmprb_clear(t);
}
