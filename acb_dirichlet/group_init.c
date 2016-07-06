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
acb_dirichlet_prime_group_init(acb_dirichlet_prime_group_struct * P, ulong p, int e)
{
    P->p = p;
    P->e = e;
    if (p == 2 || p == 4)
    {
        P->p = 2;
        nmod_init(&P->pe, 1 << e);
        if (p == 2)
        {
            P->e = 2;
            P->phi = 2;
            P->g = (1 << e) - 1;
        }
        else
        {
            P->phi = 1 << (e - 2);
            P->g = 5;
        }
    }
    else
    {
        ulong pe1;
        pe1 = n_pow(p, e - 1);
        P->phi = (p-1) * pe1;
        nmod_init(&P->pe, p * pe1);
        P->g = primitive_root_p_and_p2(p);
    }
    P->dlog = NULL;
}

static void
acb_dirichlet_group_lift_generators(acb_dirichlet_group_t G)
{
    slong k;
    acb_dirichlet_prime_group_struct * P = G->P;

    G->expo = G->phi_q = 1;
    if (G->neven)
    {
        G->phi_q = G->q_even / 2;
        G->expo = P[G->neven - 1].phi;
    }
    for (k = G->neven; k < G->num; k++)
    {
        G->phi_q *= P[k].phi;
        G->expo *= P[k].phi / n_gcd(G->expo, P[k].p - 1);
    }

    for (k = 0; k < G->num; k++)
    {
        nmod_t pe;
        ulong qpe, v;
        G->PHI[k] = G->expo / G->P[k].phi;
        /* lift generators mod q */
        /* u * p^e + v * q/p^e = 1 -> g mod q = 1 + (g-1) * v*(q/p^e) */
        pe = G->P[k].pe;
        qpe = G->q / pe.n;
        v = nmod_inv(qpe % pe.n, pe);
        /* no overflow since v * qpe < q */
        G->generators[k] = (1 + (G->P[k].g-1) * v * qpe) % G->q;
    }
}

void
acb_dirichlet_group_init(acb_dirichlet_group_t G, ulong q)
{
    slong k;
    ulong e2;
    n_factor_t fac;

    G->q = q;
    nmod_init(&G->mod, q);

    G->q_even = 1;

    for (e2 = 0; q % 2 == 0; e2++)
    {
        q /= 2;
        G->q_even *= 2;
    }

    /* warning: only factor odd part */
    n_factor_init(&fac);
    n_factor(&fac, q, 1);

    /* number of components at p=2 */
    G->neven = (e2 >= 3) ? 2 : (e2 == 2) ? 1 : 0;
    G->num = fac.num + G->neven;
    G->P = flint_malloc(G->num * sizeof(acb_dirichlet_prime_group_struct));
    G->generators = flint_malloc(G->num * sizeof(ulong));
    G->PHI = flint_malloc(G->num * sizeof(ulong));

    /* even part */
    if (G->neven >= 1)
        acb_dirichlet_prime_group_init(&G->P[0], 2, e2);
    if (G->neven == 2)
        acb_dirichlet_prime_group_init(&G->P[1], 4, e2);

    /* odd part */
    G->rad_q = 1;
    for (k = G->neven; k < G->num; k++)
    {
        ulong p, e;
        p = fac.p[k - G->neven];
        e = fac.exp[k - G->neven];
        G->rad_q *= p;
        acb_dirichlet_prime_group_init(&G->P[k], p, e);
    }
    acb_dirichlet_group_lift_generators(G);
}

void
acb_dirichlet_subgroup_init(acb_dirichlet_group_t H, const acb_dirichlet_group_t G, ulong h)
{
    slong k, j;
    int s[15]; /* selected components */

    H->q = h;
    for (k = 0, j = 0; k < G->num; k++)
    {
        ulong p = G->P[k].p;

        for(s[k] = 0; h % p == 0; h /= p)
            s[k]++;

        if (s[k] > 0)
            j++;
    }

    H->num = j;
    H->P = flint_malloc(j * sizeof(acb_dirichlet_prime_group_struct));

    j = 0;
    for (k = 0; k < G->num; k++)
    {
        if (s[k])
        {
            H->P[j] = G->P[k];
            H->P[j].e = s[k];
            j++;
        }
    }

    acb_dirichlet_group_lift_generators(H);
}
