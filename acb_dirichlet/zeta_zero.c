/*
    Copyright (C) 2010 Juan Arias de Reyna
    Copyright (C) 2019 D.H.J. Polymath

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

void
acb_dirichlet_zeta_zero(acb_t res, const fmpz_t n, slong prec)
{
    fmpz_t k;
    fmpz_init(k);
    switch (fmpz_sgn(n))
    {
        case -1:
            acb_set_d(res, 0.5);
            fmpz_neg(k, n);
            acb_dirichlet_hardy_z_zero(acb_imagref(res), k, prec);
            acb_conj(res, res);
            break;
        case 1:
            acb_set_d(res, 0.5);
            acb_dirichlet_hardy_z_zero(acb_imagref(res), n, prec);
            break;
        default:
            acb_indeterminate(res);
    }
    fmpz_clear(k);
}
