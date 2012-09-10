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

void
fmprb_div(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)
{
    long r;

    if (fmprb_contains_zero(y))
    {
        fmpr_zero(fmprb_midref(z));
        fmpr_pos_inf(fmprb_radref(z));
        return;
    }

    if (fmprb_is_exact(y))
    {
        if (fmprb_is_exact(x))
        {
            r = fmpr_div(fmprb_midref(z), fmprb_midref(x), fmprb_midref(y), prec, FMPR_RND_DOWN);
            fmpr_set_error_result(fmprb_radref(z), fmprb_midref(z), r);
        }
        else
        {
            /* (x + a) / y = x/y + a/y */

            fmpr_div(fmprb_radref(z), fmprb_radref(x), fmprb_midref(y), FMPRB_RAD_PREC, FMPR_RND_UP);
            fmpr_abs(fmprb_radref(z), fmprb_radref(z));

            r = fmpr_div(fmprb_midref(z), fmprb_midref(x), fmprb_midref(y), prec, FMPR_RND_DOWN); 
            fmpr_add_error_result(fmprb_radref(z), fmprb_radref(z),
                fmprb_midref(z), r, FMPRB_RAD_PREC, FMPR_RND_UP);
        }
    }
    else
    {
        /* x/y - (x+a)/(y-b)  =  a/(b-y) + x*b/(y*(b-y)) */
        /* where we assume y > b */

        fmpr_t t, u, by;

        fmpr_init(t);
        fmpr_init(u);
        fmpr_init(by);

        /* b - y */
        fmpr_sub(by, fmprb_radref(y), fmprb_midref(y), FMPRB_RAD_PREC, FMPR_RND_DOWN);

        /* y * (b - y) */
        fmpr_mul(t, fmprb_midref(y), by, FMPRB_RAD_PREC, FMPR_RND_DOWN);
        /* x * b */
        fmpr_mul(u, fmprb_midref(x), fmprb_radref(y), FMPRB_RAD_PREC, FMPR_RND_UP);
        /* x*b / (y*(b-y)) */
        fmpr_div(u, u, t, FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_abs(u, u);

        /* a / (b-y) */
        fmpr_div(t, fmprb_radref(x), by, FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_abs(t, t);

        fmpr_add(t, t, u, FMPRB_RAD_PREC, FMPR_RND_UP);

        r = fmpr_div(fmprb_midref(z), fmprb_midref(x), fmprb_midref(y), prec, FMPR_RND_DOWN);
        fmpr_add_error_result(fmprb_radref(z), t,
            fmprb_midref(z), r, FMPRB_RAD_PREC, FMPR_RND_UP);

        fmpr_clear(t);
        fmpr_clear(u);
        fmpr_clear(by);
    }

    fmprb_adjust(z);
}

void
fmprb_div_ui(fmprb_t z, const fmprb_t x, ulong y, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_set_ui(t, y);
    fmprb_div(z, x, t, prec);
    fmprb_clear(t);
}

void
fmprb_div_si(fmprb_t z, const fmprb_t x, long y, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_set_si(t, y);
    fmprb_div(z, x, t, prec);
    fmprb_clear(t);
}

void
fmprb_div_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_set_fmpz(t, y);
    fmprb_div(z, x, t, prec);
    fmprb_clear(t);
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
