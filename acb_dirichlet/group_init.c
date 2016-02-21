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

    Copyright (C) 2015 Jonathan Bober
    Copyright (C) 2016 Fredrik Johansson

******************************************************************************/

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
    n_factor_t fac;

    G->q = q;
    G->q_odd = q;
    G->q_even = 1;

    while (G->q_odd % 2 == 0)
    {
        G->q_odd /= 2;
        G->q_even *= 2;
    }

    n_factor_init(&fac);
    n_factor(&fac, G->q_odd, 1);

    G->num = fac.num;
    G->primes = flint_malloc(G->num * sizeof(ulong));
    G->exponents = flint_malloc(G->num * sizeof(ulong));
    G->generators = flint_malloc(G->num * sizeof(ulong));
    G->PHI = flint_malloc(G->num * sizeof(ulong));

    for (k = 0; k < G->num; k++)
    {
        G->primes[k] = fac.p[k];
        G->exponents[k] = fac.exp[k];
    }

    G->phi_q_odd = 1;
    for (k = 0; k < G->num; k++)
        G->phi_q_odd *= (G->primes[k] - 1) * n_pow(G->primes[k], G->exponents[k]-1);

    if (G->q_even == 1)
        G->phi_q = G->phi_q_odd;
    else
        G->phi_q = G->phi_q_odd * (G->q_even / 2);

    for (k = 0; k < G->num; k++)
    {
        ulong phi;

        G->generators[k] = primitive_root_p_and_p2(G->primes[k]);
        phi = n_pow(G->primes[k], G->exponents[k] - 1) * (G->primes[k] - 1);
        G->PHI[k] = G->phi_q_odd / phi;
    }
}

