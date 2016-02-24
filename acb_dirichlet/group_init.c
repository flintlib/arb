/*
    Copyright (C) 2015 Jonathan Bober
    Copyright (C) 2016 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

static ulong
primitive_root_p_and_p2(ulong p)
{
    if (p == 40487)
        return 10;

#if FLINT_BITS == 64
    if (p == UWORD(6692367337))
        return 7;

    if (p > UWORD(1000000000000))
    {
        printf("primitive root: p > 10^12 not implemented");
        abort();
    }
#endif

    return n_primitive_root_prime(p);
}

void
acb_dirichlet_group_init(acb_dirichlet_group_t G, ulong q)
{
    slong k;
    ulong e2 = 0;
    n_factor_t fac;

    G->q = q;
    G->q_odd = q;
    G->q_even = 1;

    while (G->q_odd % 2 == 0)
    {
        G->q_odd /= 2;
        G->q_even *= 2;
        e2++;
    }

    n_factor_init(&fac);
    n_factor(&fac, G->q_odd, 1);

    /* number of components at p=2 */
    G->neven = (e2 >= 3) ? 2 : (e2 == 2) ? 1 : 0;
    G->num = G->neven + fac.num;
    G->primes = flint_malloc(G->num * sizeof(ulong));
    G->exponents = flint_malloc(G->num * sizeof(ulong));
    G->generators = flint_malloc(G->num * sizeof(ulong));
    G->phi = flint_malloc(G->num * sizeof(ulong));
    G->PHI = flint_malloc(G->num * sizeof(ulong));

    /* even part */
    if (G->q_even <= 2)
    {
      G->expo = G->phi_q = 1;
    }
    else
    {
      G->phi_q = G->q_even / 2;
      G->expo = G->phi_q / 2;
    }

    for (k = 0; k < G->neven; k++)
    {
        G->primes[k] = 2;
        G->exponents[k] = (k==0) ? 1 : e2-2;
        G->generators[k] = (k==0) ? -1 : 5;
        G->phi[k] = (k==0) ? 1 : G->expo;
    }

    for (k = G->neven; k < G->num; k++)
    {
        ulong phik, p1;
        G->primes[k] = fac.p[k - G->neven];
        G->exponents[k] = fac.exp[k - G->neven];
        G->generators[k] = primitive_root_p_and_p2(G->primes[k]);
        p1 = G->primes[k] - 1;
        phik = p1 * n_pow(G->primes[k], G->exponents[k]-1);
        G->expo *= phik / n_gcd(G->expo, p1);
        G->phi_q *= phik;
        G->phi[k] = phik;
    }

    for (k = 0; k < G->num; k++)
    {
        G->PHI[k] = G->expo / G->phi[k];
        /* FIXME: generators[k] should be lifted mod q! */
    }
}
