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

    Copyright (C) 2015 Arb authors

******************************************************************************/

#include "arb.h"

void
_arb_sinc_derivative_bound(mag_t d, const arb_t x)
{
    /* |f'(x)| < min(arb_get_mag(x), 1) / 2 */
    mag_t r, one;
    mag_init(r);
    mag_init(one);
    arb_get_mag(r, x);
    mag_one(one);
    mag_min(d, r, one);
    mag_mul_2exp_si(d, d, -1);
    mag_clear(r);
    mag_clear(one);
}

void
_arb_sinc_direct(arb_t z, const arb_t x, slong prec)
{
    /* z = sin(x) / x */
    slong wp;
    arb_t y;
    wp = prec + 2;
    arb_init(y);
    arb_sin(y, x, wp);
    arb_div(z, y, x, prec);
    arb_clear(y);
}

void
arb_sinc(arb_t z, const arb_t x, slong prec)
{
    mag_t c, r;
    mag_init(c);
    mag_init(r);
    mag_set_ui_2exp_si(c, 5, -1);
    arb_get_mag_lower(r, x);
    if (mag_cmp(c, r) < 0)
    {
        _arb_sinc_direct(z, x, prec);
    }
    else if (mag_cmp_2exp_si(arb_radref(x), 1) < 0)
    {
        arb_t u;
        arb_init(u);

        /* evaluate sinc of the midpoint of x */
        if (arf_is_zero(arb_midref(x)))
        {
            arb_one(u);
        }
        else
        {
            arb_get_mid_arb(u, x);
            _arb_sinc_direct(u, u, prec);
        }

        /* account for the radius using the derivative bound */
        if (!arb_is_exact(x))
        {
            mag_t d;
            mag_init(d);
            _arb_sinc_derivative_bound(d, x);
            mag_addmul(arb_radref(u), arb_radref(x), d);
            mag_clear(d);
        }

        arb_set(z, u);
        arb_clear(u);
    }
    else
    {
        /* x has a large radius and includes points near the origin */
        arf_zero(arb_midref(z));
        mag_one(arb_radref(z));
    }

    mag_clear(c);
    mag_clear(r);
}
