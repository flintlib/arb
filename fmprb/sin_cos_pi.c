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
fmprb_sin_pi(fmprb_t y, const fmprb_t x, long prec)
{
    fmprb_t t;
    fmprb_t u;
    fmpz_t v;

    if (fmpr_cmpabs_2exp_si(fmprb_midref(t), FLINT_MAX(65536, (4*prec))) > 0)
    {
        fmpr_zero(fmprb_midref(y));
        fmpr_one(fmprb_radref(y));
        return;
    }

    fmprb_init(t);
    fmprb_init(u);
    fmpz_init(v);

    fmprb_mul_2exp_si(t, x, 1);
    fmpr_get_fmpz(v, fmprb_midref(t), FMPR_RND_NEAR);
    fmprb_sub_fmpz(t, t, v, prec);

    fmprb_const_pi(u, prec);
    fmprb_mul(t, t, u, prec);
    fmprb_mul_2exp_si(t, t, -1);

    switch (fmpz_fdiv_ui(v, 4))
    {
        case 0:
            fmprb_sin(y, t, prec);
            break;
        case 1:
            fmprb_cos(y, t, prec);
            break;
        case 2:
            fmprb_sin(y, t, prec);
            fmprb_neg(y, y);
            break;
        default:
            fmprb_cos(y, t, prec);
            fmprb_neg(y, y);
            break;
    }

    fmpz_clear(v);
    fmprb_clear(t);
    fmprb_clear(u);
}

void
fmprb_cos_pi(fmprb_t y, const fmprb_t x, long prec)
{
    fmprb_t t;
    fmprb_t u;
    fmpz_t v;

    if (fmpr_cmpabs_2exp_si(fmprb_midref(t), FLINT_MAX(65536, (4*prec))) > 0)
    {
        fmpr_zero(fmprb_midref(y));
        fmpr_one(fmprb_radref(y));
        return;
    }

    fmprb_init(t);
    fmprb_init(u);
    fmpz_init(v);

    fmprb_mul_2exp_si(t, x, 1);
    fmpr_get_fmpz(v, fmprb_midref(t), FMPR_RND_NEAR);
    fmprb_sub_fmpz(t, t, v, prec);

    fmprb_const_pi(u, prec);
    fmprb_mul(t, t, u, prec);
    fmprb_mul_2exp_si(t, t, -1);

    switch (fmpz_fdiv_ui(v, 4))
    {
        case 0:
            fmprb_cos(y, t, prec);
            break;
        case 1:
            fmprb_sin(y, t, prec);
            fmprb_neg(y, y);
            break;
        case 2:
            fmprb_cos(y, t, prec);
            fmprb_neg(y, y);
            break;
        default:
            fmprb_sin(y, t, prec);
            break;
    }

    fmpz_clear(v);
    fmprb_clear(t);
    fmprb_clear(u);
}

void
fmprb_sin_cos_pi(fmprb_t s, fmprb_t c, const fmprb_t x, long prec)
{
    fmprb_t t;
    fmprb_t u;
    fmpz_t v;

    if (fmpr_cmpabs_2exp_si(fmprb_midref(x), FLINT_MAX(65536, (4*prec))) > 0)
    {
        fmpr_zero(fmprb_midref(s));
        fmpr_one(fmprb_radref(s));
        fmpr_zero(fmprb_midref(c));
        fmpr_one(fmprb_radref(c));
        return;
    }

    fmprb_init(t);
    fmprb_init(u);
    fmpz_init(v);

    fmprb_mul_2exp_si(t, x, 1);
    fmpr_get_fmpz(v, fmprb_midref(t), FMPR_RND_NEAR);
    fmprb_sub_fmpz(t, t, v, prec);

    fmprb_const_pi(u, prec);
    fmprb_mul(t, t, u, prec);
    fmprb_mul_2exp_si(t, t, -1);

    switch (fmpz_fdiv_ui(v, 4))
    {
        case 0:
            fmprb_sin_cos(s, c, t, prec);
            break;
        case 1:
            fmprb_sin_cos(c, s, t, prec);
            fmprb_neg(c, c);
            break;
        case 2:
            fmprb_sin_cos(s, c, t, prec);
            fmprb_neg(s, s);
            fmprb_neg(c, c);
            break;
        default:
            fmprb_sin_cos(c, s, t, prec);
            fmprb_neg(s, s);
            break;
    }

    fmpz_clear(v);
    fmprb_clear(t);
    fmprb_clear(u);
}

