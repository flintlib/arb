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

long
_fmpr_sin(fmpr_t y, const fmpr_t x, long prec, fmpr_rnd_t rnd)
{
    if (fmpr_is_special(x))
    {
        if (fmpr_is_zero(x))
            fmpr_zero(y);
        else
            fmpr_nan(y);

        return FMPR_RESULT_EXACT;
    }
    else
    {
        long r;
        CALL_MPFR_FUNC(r, mpfr_sin, y, x, prec, rnd);
        return r;
    }
}

long
_fmpr_cos(fmpr_t y, const fmpr_t x, long prec, fmpr_rnd_t rnd)
{
    if (fmpr_is_special(x))
    {
        if (fmpr_is_zero(x))
            fmpr_one(y);
        else
            fmpr_nan(y);

        return FMPR_RESULT_EXACT;
    }
    else
    {
        long r;
        CALL_MPFR_FUNC(r, mpfr_cos, y, x, prec, rnd);
        return r;
    }
}

void
fmprb_sin(fmprb_t z, const fmprb_t x, long prec)
{
    long r;

    if (fmprb_is_exact(x))
    {
        r = _fmpr_sin(fmprb_midref(z), fmprb_midref(x), prec, FMPR_RND_DOWN);
        fmpr_set_error_result(fmprb_radref(z), fmprb_midref(z), r);
    }
    else
    {
        fmpr_t t;
        fmpr_init(t);
        fmpr_set(t, fmprb_radref(x));

        r = _fmpr_sin(fmprb_midref(z), fmprb_midref(x), prec, FMPR_RND_DOWN);
        fmpr_add_error_result(fmprb_radref(z), t, fmprb_midref(z), r,
            FMPRB_RAD_PREC, FMPR_RND_UP);

        fmpr_clear(t);
    }

    fmprb_adjust(z);
}

void
fmprb_cos(fmprb_t z, const fmprb_t x, long prec)
{
    long r;

    if (fmprb_is_exact(x))
    {
        r = _fmpr_cos(fmprb_midref(z), fmprb_midref(x), prec, FMPR_RND_DOWN);
        fmpr_set_error_result(fmprb_radref(z), fmprb_midref(z), r);
    }
    else
    {
        fmpr_t t;
        fmpr_init(t);
        fmpr_set(t, fmprb_radref(x));

        r = _fmpr_cos(fmprb_midref(z), fmprb_midref(x), prec, FMPR_RND_DOWN);
        fmpr_add_error_result(fmprb_radref(z), t, fmprb_midref(z), r,
            FMPRB_RAD_PREC, FMPR_RND_UP);

        fmpr_clear(t);
    }

    fmprb_adjust(z);
}

void
fmprb_sin_cos(fmprb_t s, fmprb_t c, const fmprb_t x, long prec)
{
    fmprb_sin(s, x, prec);
    fmprb_cos(c, x, prec);
}
