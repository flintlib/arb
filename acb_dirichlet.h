/*
    Copyright (C) 2015 Jonathan Bober
    Copyright (C) 2016 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef ACB_DIRICHLET_H
#define ACB_DIRICHLET_H

#include "acb.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    ulong q;                /* modulus */
    ulong q_even;           /* even part of modulus */
    ulong q_odd;            /* odd part of modulus */
    ulong phi_q;            /* phi(q) */
    ulong phi_q_odd;        /* phi(q_odd) */
    slong num;              /* number of odd prime factors */
    ulong * primes;         /* odd primes p[k] */
    ulong * exponents;      /* exponents e[k] */
    ulong * generators;     /* generator for each odd prime p[k] */
    ulong * PHI;            /* PHI(k) = phi(q_odd) / phi(p[k]^e[k])  */
}
acb_dirichlet_group_struct;

typedef acb_dirichlet_group_struct acb_dirichlet_group_t[1];

void _acb_dirichlet_euler_product_real_ui(arb_t res, ulong s,
    const signed char * chi, int mod, int reciprocal, slong prec);

void acb_dirichlet_eta(acb_t res, const acb_t s, slong prec);

void acb_dirichlet_group_init(acb_dirichlet_group_t G, ulong q);

void acb_dirichlet_group_clear(acb_dirichlet_group_t G);

void acb_dirichlet_chi(acb_t res, const acb_dirichlet_group_t G, ulong m, ulong n, slong prec);

#ifdef __cplusplus
}
#endif

#endif

