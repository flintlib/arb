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
acb_dirichlet_arb_theta_naive(acb_t res, const arb_t x, int parity, const ulong *a, const acb_dirichlet_powers_t z, slong len, slong prec)
{
    slong k;
    arb_t xk2, dx, x2;
    acb_t zk;
    arb_init(xk2);
    arb_init(dx);
    arb_init(x2);
    acb_init(zk);

    arb_set(dx, x);
    arb_set(xk2, dx);
    arb_mul(x2, dx, dx, prec);

    acb_set_arb(res, xk2);

    /* TODO: reduce prec */
    for (k = 2; k < len; k++)
    {
        arb_mul(dx, dx, x2, prec);
        arb_mul(xk2, xk2, dx, prec);
        if (a[k] != ACB_DIRICHLET_CHI_NULL)
        {
            acb_dirichlet_power(zk, z, a[k], prec);
            if (parity)
                acb_mul_si(zk, zk, k, prec);
            acb_addmul_arb(res, zk, xk2, prec);
        }
    }

    arb_clear(xk2);
    arb_clear(x2);
    arb_clear(dx);
    acb_clear(zk);
}
