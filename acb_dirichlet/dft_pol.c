/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_dirichlet.h"

/* all roots are already computed */
void
_acb_dirichlet_dft_pol(acb_ptr w, acb_srcptr v, slong dv, acb_srcptr z, slong dz, slong len, slong prec)
{
    slong i, j;
    acb_ptr wi;
    acb_srcptr vj;

    for (i = 0, wi = w; i < len; i++, wi++)
    {
        acb_zero(wi);
        for (j = 0, vj = v; j < len; j++, vj += dv)
            acb_addmul(wi, vj, z + dz * (i * j % len), prec);
    }
}

void
acb_dirichlet_dft_pol(acb_ptr w, acb_srcptr v, slong len, slong prec)
{
    acb_ptr z;
    z = _acb_vec_init(len);
    acb_dirichlet_vec_nth_roots(z, len, prec);
    _acb_dirichlet_dft_pol(w, v, 1, z, 1, len, prec);
    _acb_vec_clear(z, len);
}
