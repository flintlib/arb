/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"
#include "acb_poly.h"

void
acb_dirichlet_power(acb_t z, const acb_dirichlet_powers_t t, ulong n, slong prec)
{
    if (!t->depth)
    {
        acb_pow_ui(z, t->z, n, prec);
    }
    else
    {
        slong k;
        ulong r;
        r = n % t->size;
        n = n / t->size;
        acb_set(z, t->Z[0] + r);
        for (k = 1; k < t->depth && n; k++)
        {
            r = n % t->size;
            n = n / t->size;
            acb_mul(z, z, t->Z[k] + r, prec);
        }
    }
}
