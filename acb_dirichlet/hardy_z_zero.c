/*
    Copyright (C) 2019 D.H.J. Polymath

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"
#include "arb_calc.h"

static void
_acb_set_arf(acb_t res, const arf_t t)
{
    acb_zero(res);
    arb_set_arf(acb_realref(res), t);
}

int
_acb_dirichlet_definite_hardy_z(arb_t res, const arf_t t, slong *pprec)
{
    int msign;
    acb_t z;
    acb_init(z);
    while (1)
    {
        _acb_set_arf(z, t);
        acb_dirichlet_hardy_z(z, z, NULL, NULL, 1, *pprec);
        msign = arb_sgn_nonzero(acb_realref(z));
        if (msign)
        {
            break;
        }
        *pprec *= 2;
    }
    acb_get_real(res, z);
    acb_clear(z);
    return msign;
}

static int
_partition_hardy_z(arf_interval_t L, arf_interval_t R,
        const arf_interval_t block, slong *pprec)
{
    arb_t t;
    arf_t u;
    int msign;

    arb_init(t);
    arf_init(u);

    /* Compute the midpoint */
    arf_add(u, &block->a, &block->b, ARF_PREC_EXACT, ARF_RND_DOWN);
    arf_mul_2exp_si(u, u, -1);

    /* Evaluate and get sign at midpoint */
    msign = 0;
    msign = _acb_dirichlet_definite_hardy_z(t, u, pprec);

    /* L, R = block, split at midpoint */
    arf_set(&L->a, &block->a);
    arf_set(&R->b, &block->b);
    arf_set(&L->b, u);
    arf_set(&R->a, u);

    arb_clear(t);
    arf_clear(u);

    return msign;
}

static void
_refine_hardy_z_zero_bisect(arf_interval_t res,
        const arf_interval_t start, slong iter)
{
    int asign, bsign, msign;
    slong i, prec = 8;
    arf_interval_t t, u;
    arb_t m, v;

    arf_interval_init(t);
    arf_interval_init(u);
    arb_init(m);
    arb_init(v);

    asign = _acb_dirichlet_definite_hardy_z(v, &start->a, &prec);
    bsign = _acb_dirichlet_definite_hardy_z(v, &start->b, &prec);
    if (asign == bsign)
    {
        flint_printf("isolate a zero before bisecting the interval\n");
        flint_abort();
    }
    arf_interval_set(res, start);
    for (i = 0; i < iter; i++)
    {
        msign = _partition_hardy_z(t, u, res, &prec);
        if (msign == asign)
            arf_interval_swap(res, u);
        else
            arf_interval_swap(res, t);
    }

    arf_interval_clear(t);
    arf_interval_clear(u);
    arb_clear(m);
    arb_clear(v);
}

void
acb_dirichlet_hardy_z_zero(arb_t res, const fmpz_t n, slong prec)
{
    arf_interval_t r, s;
    mag_t m;
    slong bits;

    if (fmpz_cmp_si(n, 1) < 0)
    {
        flint_printf("n must be positive\n");
        flint_abort();
    }

    arf_interval_init(r);
    arf_interval_init(s);
    mag_init(m);

    acb_dirichlet_isolate_hardy_z_zero(&r->a, &r->b, n);

    arf_get_mag(m, &r->b);
    bits = FLINT_MAX(0, mag_get_d_log2_approx(m)) + 8;
    arb_set_interval_arf(res, &r->a, &r->b, bits);
    bits = arb_rel_accuracy_bits(res);

    if (bits < prec)
    {
        _refine_hardy_z_zero_bisect(s, r, prec - bits);
        arb_set_interval_arf(res, &s->a, &s->b, prec);
    }

    arb_set_round(res, res, prec);

    arf_interval_clear(r);
    arf_interval_clear(s);
    mag_clear(m);
}
