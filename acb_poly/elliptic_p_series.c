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

#include "acb_poly.h"
#include "acb_modular.h"

void
_acb_poly_elliptic_p_series(acb_ptr res, acb_srcptr z, long zlen, const acb_t tau, long len, long prec)
{
    acb_ptr t, u;

    zlen = FLINT_MIN(zlen, len);

    t = _acb_vec_init(len);
    u = _acb_vec_init(len);

    acb_modular_elliptic_p_zpx(t, z, tau, len, prec);

    /* compose with nonconstant part */
    acb_zero(u);
    _acb_vec_set(u + 1, z + 1, zlen - 1);
    _acb_poly_compose_series(res, t, len, u, zlen, len, prec);

    _acb_vec_clear(t, len);
    _acb_vec_clear(u, len);
}

void
acb_poly_elliptic_p_series(acb_poly_t res, const acb_poly_t z, const acb_t tau, long n, long prec)
{
    if (n == 0)
    {
        acb_poly_zero(res);
        return;
    }

    acb_poly_fit_length(res, n);

    if (z->length == 0)
    {
        acb_t t;
        acb_init(t);
        _acb_poly_elliptic_p_series(res->coeffs, t, 1, tau, n, prec);
        acb_clear(t);
    }
    else
    {
        _acb_poly_elliptic_p_series(res->coeffs, z->coeffs, z->length, tau, n, prec);
    }

    _acb_poly_set_length(res, n);
    _acb_poly_normalise(res);
}
