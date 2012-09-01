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
fmprb_mul(fmprb_t z, const fmprb_t x, const fmprb_t y, long prec)
{
    long r;

    if (fmprb_is_exact(x))
    {
        if (fmprb_is_exact(y))
        {
            r = fmpr_mul(fmprb_midref(z), fmprb_midref(x), fmprb_midref(y), prec, FMPR_RND_DOWN);
            fmpr_set_error_result(fmprb_radref(z), fmprb_midref(z), r);
        }
        else
        {
            /* x * (y+b) = x*y + x*b */
            fmpr_mul(fmprb_radref(z), fmprb_midref(x), fmprb_radref(y), FMPRB_RAD_PREC, FMPR_RND_UP);
            fmpr_abs(fmprb_radref(z), fmprb_radref(z));

            r = fmpr_mul(fmprb_midref(z), fmprb_midref(x), fmprb_midref(y), prec, FMPR_RND_DOWN); 
            fmpr_add_error_result(fmprb_radref(z), fmprb_radref(z),
                fmprb_midref(z), r, FMPRB_RAD_PREC, FMPR_RND_UP);
        }
    }
    else
    {
        if (fmprb_is_exact(x))
        {
            /* (x+a) * y = x*y + y*a */
            fmpr_mul(fmprb_radref(z), fmprb_midref(y), fmprb_radref(x), FMPRB_RAD_PREC, FMPR_RND_UP);
            fmpr_abs(fmprb_radref(z), fmprb_radref(z));

            r = fmpr_mul(fmprb_midref(z), fmprb_midref(x), fmprb_midref(y), prec, FMPR_RND_DOWN);
            fmpr_add_error_result(fmprb_radref(z), fmprb_radref(z),
                fmprb_midref(z), r, FMPRB_RAD_PREC, FMPR_RND_UP);
        }
        else
        {
            /* (x+a)*(y+b) = x*y + x*b + y*a + a*b*/
            fmpr_t t, u;
            fmpr_init(t);
            fmpr_init(u);

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

    fmprb_adjust(z);
}

void
fmprb_mul_ui(fmprb_t z, const fmprb_t x, ulong y, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_set_ui(t, y);
    fmprb_mul(z, x, t, prec);
    fmprb_clear(t);
}

void
fmprb_mul_si(fmprb_t z, const fmprb_t x, long y, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_set_si(t, y);
    fmprb_mul(z, x, t, prec);
    fmprb_clear(t);
}

void
fmprb_mul_fmpz(fmprb_t z, const fmprb_t x, const fmpz_t y, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_set_fmpz(t, y);
    fmprb_mul(z, x, t, prec);
    fmprb_clear(t);
}
