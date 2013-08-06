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
fmprb_mul_fmpr_naive(fmprb_t z, const fmprb_t x, const fmpr_t y, long prec)
{
    /* (x+a) * y = x*y + y*a */
    if (fmpr_is_pos_inf(fmprb_radref(x)))
    {
        fmpr_mul(fmprb_midref(z), fmprb_midref(x), y, prec, FMPR_RND_DOWN);
        fmpr_pos_inf(fmprb_radref(z));
    }
    else if (fmpr_is_zero(fmprb_radref(x)))
    {
        long r;

        r = fmpr_mul(fmprb_midref(z), fmprb_midref(x), y, prec, FMPR_RND_DOWN);
        fmpr_set_error_result(fmprb_radref(z), fmprb_midref(z), r);
    }
    else
    {
        long r;

        fmpr_mul(fmprb_radref(z), fmprb_radref(x), y, FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_abs(fmprb_radref(z), fmprb_radref(z));

        r = fmpr_mul(fmprb_midref(z), fmprb_midref(x), y, prec, FMPR_RND_DOWN);
        fmpr_add_error_result(fmprb_radref(z), fmprb_radref(z),
            fmprb_midref(z), r, FMPRB_RAD_PREC, FMPR_RND_UP);
    }
}

void
fmprb_mul_main_naive(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)
{
    if (fmpr_is_pos_inf(fmprb_radref(x)) || fmpr_is_pos_inf(fmprb_radref(y)))
    {
        fmpr_mul(fmprb_midref(z), fmprb_midref(x), fmprb_midref(y), prec, FMPR_RND_DOWN);
        fmpr_pos_inf(fmprb_radref(z));
    }
    else
    {
        fmpr_t t, u;
        long r;

        fmpr_init(t);
        fmpr_init(u);

        /* (x+a)*(y+b) = x*y + x*b + y*a + a*b*/

        fmpr_mul(t, fmprb_midref(x), fmprb_radref(y), FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_abs(t, t);

        fmpr_mul(u, fmprb_midref(y), fmprb_radref(x), FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_abs(u, u);

        fmpr_add(t, t, u, FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_addmul(t, fmprb_radref(x), fmprb_radref(y), FMPRB_RAD_PREC, FMPR_RND_UP);

        r = fmpr_mul(fmprb_midref(z), fmprb_midref(x), fmprb_midref(y), prec, FMPR_RND_DOWN);
        fmpr_add_error_result(fmprb_radref(z), t,
            fmprb_midref(z), r, FMPRB_RAD_PREC, FMPR_RND_UP);

        fmpr_clear(t);
        fmpr_clear(u);
    }
}

void
fmprb_mul_naive(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)
{
    if (fmprb_is_exact(x))
        fmprb_mul_fmpr_naive(z, y, fmprb_midref(x), prec);
    else if (fmprb_is_exact(y))
        fmprb_mul_fmpr_naive(z, x, fmprb_midref(y), prec);
    else
        fmprb_mul_main_naive(z, x, y, prec);

    fmprb_adjust(z);
}

