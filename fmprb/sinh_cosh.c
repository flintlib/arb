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
_fmpr_sinh(fmpr_t y, const fmpr_t x, long prec, fmpr_rnd_t rnd)
{
    if (fmpr_is_special(x))
    {
        if (fmpr_is_zero(x))
            fmpr_zero(y);
        else if (fmpr_is_inf(x))
            fmpr_set(y, x);
        else
            fmpr_nan(y);

        return FMPR_RESULT_EXACT;
    }
    else
    {
        long r;
        CALL_MPFR_FUNC(r, mpfr_sinh, y, x, prec, rnd);
        return r;
    }
}

long
_fmpr_cosh(fmpr_t y, const fmpr_t x, long prec, fmpr_rnd_t rnd)
{
    if (fmpr_is_special(x))
    {
        if (fmpr_is_zero(x))
            fmpr_one(y);
        else if (fmpr_is_inf(x))
            fmpr_abs(y, x);
        else
            fmpr_nan(y);

        return FMPR_RESULT_EXACT;
    }
    else
    {
        long r;
        CALL_MPFR_FUNC(r, mpfr_cosh, y, x, prec, rnd);
        return r;
    }
}

void
_fmpr_sinh_cosh(long * r1, long * r2, fmpr_t s, fmpr_t c, const fmpr_t x, long prec, fmpr_rnd_t rnd)
{
    if (fmpr_is_special(x))
    {
        if (fmpr_is_zero(x))
        {
            fmpr_zero(s);
            fmpr_one(c);
        }
        else if (fmpr_is_pos_inf(x))
        {
            fmpr_pos_inf(s);
            fmpr_pos_inf(c);
        }
        else if (fmpr_is_neg_inf(x))
        {
            fmpr_neg_inf(s);
            fmpr_pos_inf(c);
        }
        else
        {
            fmpr_nan(s);
            fmpr_nan(c);
        }

        *r1 = *r2 = FMPR_RESULT_EXACT;
    }
    else
    {
        CALL_MPFR_FUNC_2X1(*r1, *r2, mpfr_sinh_cosh, s, c, x, prec, rnd);
    }
}

void
fmprb_sinh(fmprb_t z, const fmprb_t x, long prec)
{
    long r;

    if (fmprb_is_exact(x))
    {
        r = _fmpr_sinh(fmprb_midref(z), fmprb_midref(x), prec, FMPR_RND_DOWN);
        fmpr_set_error_result(fmprb_radref(z), fmprb_midref(z), r);
    }
    else
    {
        fmpr_t t;
        fmpr_init(t);

        /* the propagated error is bounded by sup(cosh([a,b]))*r */
        if (fmpr_sgn(fmprb_midref(x)) >= 0)
            fmpr_add(t, fmprb_midref(x), fmprb_radref(x), FMPRB_RAD_PREC, FMPR_RND_UP);
        else
            fmpr_sub(t, fmprb_radref(x), fmprb_midref(x), FMPRB_RAD_PREC, FMPR_RND_UP);

        _fmpr_cosh(t, t, FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_mul(t, t, fmprb_radref(x), FMPRB_RAD_PREC, FMPR_RND_UP);

        r = _fmpr_sinh(fmprb_midref(z), fmprb_midref(x), prec, FMPR_RND_DOWN);
        fmpr_add_error_result(fmprb_radref(z), t, fmprb_midref(z), r,
            FMPRB_RAD_PREC, FMPR_RND_UP);

        fmpr_clear(t);
    }

    fmprb_adjust(z);
}

void
fmprb_cosh(fmprb_t z, const fmprb_t x, long prec)
{
    long r;

    if (fmprb_is_exact(x))
    {
        r = _fmpr_cosh(fmprb_midref(z), fmprb_midref(x), prec, FMPR_RND_DOWN);
        fmpr_set_error_result(fmprb_radref(z), fmprb_midref(z), r);
    }
    else
    {
        fmpr_t t;
        fmpr_init(t);

        /* the propagated error is bounded by sup(|sinh([a,b])|)*r */
        if (fmpr_sgn(fmprb_midref(x)) >= 0)
            fmpr_add(t, fmprb_midref(x), fmprb_radref(x), FMPRB_RAD_PREC, FMPR_RND_UP);
        else
            fmpr_sub(t, fmprb_radref(x), fmprb_midref(x), FMPRB_RAD_PREC, FMPR_RND_UP);

        _fmpr_sinh(t, t, FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_mul(t, t, fmprb_radref(x), FMPRB_RAD_PREC, FMPR_RND_UP);

        r = _fmpr_cosh(fmprb_midref(z), fmprb_midref(x), prec, FMPR_RND_DOWN);
        fmpr_add_error_result(fmprb_radref(z), t, fmprb_midref(z), r,
            FMPRB_RAD_PREC, FMPR_RND_UP);

        fmpr_clear(t);
    }

    fmprb_adjust(z);
}

void
fmprb_sinh_cosh(fmprb_t s, fmprb_t c, const fmprb_t x, long prec)
{
    long r1, r2;

    if (fmprb_is_exact(x))
    {
        _fmpr_sinh_cosh(&r1, &r2, fmprb_midref(s), fmprb_midref(c), fmprb_midref(x), prec, FMPR_RND_DOWN);
        fmpr_set_error_result(fmprb_radref(s), fmprb_midref(s), r1);
        fmpr_set_error_result(fmprb_radref(c), fmprb_midref(c), r2);
    }
    else
    {
        fmpr_t t, u;
        fmpr_init(t);
        fmpr_init(u);

        if (fmpr_sgn(fmprb_midref(x)) >= 0)
            fmpr_add(t, fmprb_midref(x), fmprb_radref(x), FMPRB_RAD_PREC, FMPR_RND_UP);
        else
            fmpr_sub(t, fmprb_radref(x), fmprb_midref(x), FMPRB_RAD_PREC, FMPR_RND_UP);

        _fmpr_sinh_cosh(&r1, &r2, t, u, t, FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_mul(t, t, fmprb_radref(x), FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_mul(u, u, fmprb_radref(x), FMPRB_RAD_PREC, FMPR_RND_UP);

        _fmpr_sinh_cosh(&r1, &r2, fmprb_midref(s), fmprb_midref(c), fmprb_midref(x), prec, FMPR_RND_DOWN);

        fmpr_add_error_result(fmprb_radref(s), u, fmprb_midref(s), r1, FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_add_error_result(fmprb_radref(c), t, fmprb_midref(c), r2, FMPRB_RAD_PREC, FMPR_RND_UP);

        fmpr_clear(t);
        fmpr_clear(u);
    }

    fmprb_adjust(s);
    fmprb_adjust(c);
}

