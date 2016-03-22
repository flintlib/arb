/*
    Copyright (C) 2015 Jonathan Bober
    Copyright (C) 2016 Fredrik Johansson
    Copyright (C) 2016 Pascal Molin

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
    G->primepowers = flint_malloc(G->num * sizeof(ulong));
    G->generators = flint_malloc(G->num * sizeof(ulong));
    G->phi = flint_malloc(G->num * sizeof(ulong));
    G->PHI = flint_malloc(G->num * sizeof(ulong));

    /* even part */
    G->expo = G->phi_q = 1;
    if (G->neven >= 1)
    {
        G->primes[0] = 2;
        G->exponents[0] = 2;
        G->phi[0] = 2;
        G->primepowers[0] = G->q_even;
        G->generators[0] = G->q_even-1;
        G->expo = 2;
        G->phi_q = 2;
    }
    if (G->neven == 2)
    {
        G->primes[1] = 2;
        G->exponents[1] = e2;
        G->phi[1] = G->q_even / 4;
        G->primepowers[1] = G->q_even;
        G->generators[1] = 5;
        G->expo = G->phi[1];
        G->phi_q = G->q_even / 2;
    }
    /* odd part */
    for (k = G->neven; k < G->num; k++)
    {
        ulong p1, pe1;
        G->primes[k] = fac.p[k - G->neven];
        G->exponents[k] = fac.exp[k - G->neven];
        p1 = G->primes[k] - 1;
        pe1 = n_pow(G->primes[k], G->exponents[k]-1);
        G->phi[k] = p1 * pe1;
        G->primepowers[k] = pe1 * G->primes[k];
        G->generators[k] = primitive_root_p_and_p2(G->primes[k]);
        G->expo *= G->phi[k] / n_gcd(G->expo, p1);
        G->phi_q *= G->phi[k];
    }
    /* generic odd+even */
    for (k = 0; k < G->num; k++)
    {
        ulong pe, qpe, v;
        G->PHI[k] = G->expo / G->phi[k];
        /* lift generators mod q */
        /* u * p^e + v * q/p^e = 1 -> g mod q = 1 + (g-1) * v*(q/p^e) */
        pe = G->primepowers[k];
        qpe = q / pe;
        v = n_invmod(qpe % pe, pe);
        /* no overflow since v * qpe < q */
        G->generators[k] = (1 + (G->generators[k]-1) * v * qpe) % q;
    }
}
