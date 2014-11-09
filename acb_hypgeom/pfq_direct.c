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

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include "acb_hypgeom.h"

void
acb_hypgeom_pfq_direct(acb_t res, acb_srcptr a, long p, acb_srcptr b, long q,
    const acb_t z, long n, long prec)
{
    acb_t s, t;
    mag_t err, C;

    acb_init(s);
    acb_init(t);
    mag_init(err);
    mag_init(C);

    if (n < 0)
        n = acb_hypgeom_pfq_choose_n(a, p, b, q, z, prec);

    acb_hypgeom_pfq_sum(s, t, a, p, b, q, z, n, prec);

    if (!acb_is_zero(t))
    {
        acb_hypgeom_pfq_bound_factor(C, a, p, b, q, z, n);
        acb_get_mag(err, t);
        mag_mul(err, err, C);

        if (_acb_vec_is_real(a, p) && _acb_vec_is_real(b, q) && acb_is_real(z))
            arb_add_error_mag(acb_realref(s), err);
        else
            acb_add_error_mag(s, err);
    }

    acb_set(res, s);

    acb_clear(s);
    acb_clear(t);
    mag_clear(err);
    mag_clear(C);
}

