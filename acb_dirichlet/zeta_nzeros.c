/*
    Copyright (C) 2019 D.H.J Polymath

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

static slong
_nzeros_small(const arf_t t)
{
    slong n = 14;
    slong interlace[14] = {
        14, 18, 23, 28, 32, 35, 38, 42, 45, 49, 52, 54, 58, 60};
    if (arf_cmp_si(t, interlace[n-1]) >= 0)
    {
        flint_printf("error: requires t < %ld\n", interlace[n-1]);
        flint_abort();
    }
    else if (arf_cmp_si(t, interlace[0]) < 0)
    {
        return 0;
    }
    else
    {
        acb_t z;
        slong k, sa, sb;
        slong prec = 4;
        acb_init(z);
        while (1)
        {
            acb_zero(z);
            arb_set_arf(acb_realref(z), t);
            acb_dirichlet_hardy_z(z, z, NULL, NULL, 1, prec);
            if (!acb_contains_zero(z))
            {
                break;
            }
            prec *= 2;
        }
        sa = arb_is_positive(acb_realref(z)) ? 1 : -1;
        acb_clear(z);
        for (k = 1; k < n; k++)
        {
            if (arf_cmp_si(t, interlace[k]) < 0)
            {
                sb = (k % 2) ? 1 : -1;
                return (sa == sb) ? k : k-1;
            }
        }
        flint_printf("error: _nzeros_small internal error\n");
        flint_abort();
    }
    return -1;
}

void
acb_dirichlet_zeta_nzeros(arb_t res, const arb_t t, slong prec)
{
    arb_t c;
    arb_init(c);
    arb_set_si(c, 60);
    if (arb_lt(t, c))
    {
        if (arb_is_exact(t))
        {
            arb_set_si(res, _nzeros_small(arb_midref(t)));
        }
        else
        {
            arf_t lb, ub;
            slong a, b;
            arf_init(lb);
            arf_init(ub);
            arb_get_lbound_arf(lb, t, prec);
            arb_get_ubound_arf(ub, t, prec);
            a = _nzeros_small(lb);
            b = _nzeros_small(ub);
            arb_set_ui(res, a+b);
            mag_set_ui(arb_radref(res), b-a);
            arb_mul_2exp_si(res, res, -1);
            arf_clear(lb);
            arf_clear(ub);
        }
    }
    else
    {
        /* N(t) = theta(t)/pi + 1 + S(t) */
        arb_t s, pi;
        acb_t z;
        fmpz_t n;
        arb_init(s);
        arb_init(pi);
        acb_init(z);
        fmpz_init(n);
        arb_const_pi(pi, prec);
        acb_dirichlet_backlund_s(s, t, prec);
        acb_set_arb(z, t);
        acb_dirichlet_hardy_theta(z, z, NULL, NULL, 1, prec);
        arb_div(res, acb_realref(z), pi, prec);
        arb_add(res, res, s, prec);
        arb_add_ui(res, res, 1, prec);
        arb_clear(s);
        arb_clear(pi);
        acb_clear(z);
        fmpz_clear(n);
    }
    arb_clear(c);
}
