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

    Copyright (C) 2016 Pascal Molin

******************************************************************************/

#include "acb_poly.h"
#include "acb_dirichlet.h"

acb_ptr
acb_roots_init(slong len, slong prec)
{
    acb_t zeta;
    acb_ptr z;
    acb_init(zeta);
    acb_dirichlet_nth_root(zeta, len, prec);
    z = _acb_vec_init(len);
    _acb_vec_set_powers(z, zeta, len, prec);
    acb_clear(zeta);
    return z;
}

void
_acb_dirichlet_dft_pol(acb_ptr w, acb_srcptr v, acb_srcptr z, slong len, slong prec)
{
    _acb_poly_evaluate_vec_fast(w, v, len, z, len, prec);
}

void
acb_dirichlet_dft_pol(acb_ptr w, acb_srcptr v, slong len, slong prec)
{
    acb_ptr z;
    z = acb_roots_init(len, prec);
    _acb_dirichlet_dft_pol(w, v, z, len, prec);
    _acb_vec_clear(z, len);
}
