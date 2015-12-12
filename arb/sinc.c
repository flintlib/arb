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
#include "acb.h"
#include "acb_hypgeom.h"

void
_arb_sinc_0f1(arb_t res, const arb_t x, slong prec)
{
    acb_struct bb[2];
    arb_t y;
    acb_t z, w;

    arb_init(y);
    acb_init(bb + 0);
    acb_init(bb + 1);
    acb_init(z);
    acb_init(w);

    /* a = 3/2 */
    acb_set_ui(bb + 0, 3);
    acb_mul_2exp_si(bb + 0, bb + 0, -1);

    /* z = -(x/2)^2 */
    arb_mul_2exp_si(y, x, -1);
    arb_mul(y, y, y, prec);
    arb_neg(y, y);
    acb_set_arb(z, y);

    /* res = 0F1(a, z) */
    acb_one(bb + 1);
    acb_hypgeom_pfq_direct(w, NULL, 0, bb, 2, z, -1, prec);
    arb_set(res, acb_realref(w));

    arb_clear(y);
    acb_clear(bb+0);
    acb_clear(bb+1);
    acb_clear(z);
    acb_clear(w);
}

void
_arb_sinc_dbound(arb_t res, const arb_t x, slong prec)
{
    mag_t r;
    arb_t a;

    mag_init(r);
    arb_init(a);

    /* compute sinc of the midpoint, with error due to inexact sin and div */
    if (arf_is_zero(arb_midref(x)))
    {
        arb_one(a);
    }
    else
    {
        arb_set_arf(a, arb_midref(x));
        arb_sin(res, a, prec);
        arb_div(res, res, a, prec);
    }

    /* add radius error using a global bound of the derivative of sinc */
    /* |sinc(x)'| < 1/2 */
    mag_mul_2exp_si(r, arb_radref(x), -1);
    mag_add(arb_radref(res), arb_radref(res), r);

    mag_clear(r);
    arb_clear(a);
}

int
_mag_cmp_ui(const mag_t x, ulong y)
{
    int res;
    mag_t z;
    mag_init(z);
    mag_set_ui(z, y);
    res = mag_cmp(x, z);
    mag_clear(z);
    return res;
}

void
_mag_addmul_ui(mag_t z, const mag_t x, ulong y)
{
    mag_t w;
    mag_init(w);
    mag_mul_ui(w, x, y);
    mag_add(z, z, w);
    mag_clear(w);
}

int
_arb_sinc_use_0f1(const arb_t x)
{
    /* 8|m| + 5r < 10 */
    int res;

    mag_t lhs;
    mag_init(lhs);
    arf_get_mag(lhs, arb_midref(x));
    mag_mul_ui(lhs, lhs, 8);
    _mag_addmul_ui(lhs, arb_radref(x), 5);

    res = (_mag_cmp_ui(lhs, 10) < 0);

    mag_clear(lhs);
    return res;
}

void
arb_sinc(arb_t res, const arb_t x, slong prec)
{
    arb_t b;
    arb_init(b);
    mag_set_ui_2exp_si(arb_radref(b), 5, -1);
    if (arb_overlaps(b, x))
    {
        if (_arb_sinc_use_0f1(x))
        {
            _arb_sinc_0f1(res, x, prec);
        }
        else if (mag_cmp_2exp_si(arb_radref(x), 1) < 0)
        {
            _arb_sinc_dbound(res, x, prec);
        }
        else
        {
            arf_zero(arb_midref(res));
            mag_one(arb_radref(res));
        }
    }
    else
    {
        arb_sin(res, x, prec);
        arb_div(res, res, x, prec);
    }
    arb_clear(b);
}
