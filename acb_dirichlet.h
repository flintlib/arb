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
    ulong phi_q;            /* phi(q) = group size */
    ulong expo;             /* group exponent = lcm(phi(q_even), phi(p[k]^e[k]) ) */
    slong neven;            /* number of even components (in 0,1,2)*/
    slong num;              /* number of prime components (even + odd) */
    ulong * primes;         /* primes p[k] */
    ulong * exponents;      /* exponents e[k] */
    ulong * primepowers;    /* powers p[k]^[k] */
    ulong * generators;     /* generator for each prime p[k] lifted mod q */
    ulong * phi;            /* phi(k) = phi(p[k]^e[k])       */
    ulong * PHI;            /* PHI(k) = expo / phi(k)        */
}
acb_dirichlet_group_struct;

typedef acb_dirichlet_group_struct acb_dirichlet_group_t[1];

/* elements of the group, keep both number and log */
typedef struct
{
  ulong n;           /* number */
  ulong * log;       /* s.t. prod generators[k]^log[k] = number */
}
acb_conrey_struct;

typedef acb_conrey_struct acb_conrey_t[1];

void _acb_dirichlet_euler_product_real_ui(arb_t res, ulong s,
    const signed char * chi, int mod, int reciprocal, slong prec);

void acb_dirichlet_eta(acb_t res, const acb_t s, slong prec);

void acb_dirichlet_group_init(acb_dirichlet_group_t G, ulong q);

void acb_dirichlet_group_clear(acb_dirichlet_group_t G);

void acb_conrey_init(acb_conrey_t x, const acb_dirichlet_group_t G);

void acb_conrey_clear(acb_conrey_t x);

void acb_conrey_one(acb_conrey_t x, const acb_dirichlet_group_t G);

void acb_conrey_first_primitive(acb_conrey_t x, const acb_dirichlet_group_t G);

void acb_conrey_log(acb_conrey_t x, const acb_dirichlet_group_t G, ulong m);

int acb_conrey_next(acb_conrey_t x, const acb_dirichlet_group_t G);

int acb_conrey_next_primitive(acb_conrey_t x, const acb_dirichlet_group_t G);

long n_dirichlet_chi_conrey(const acb_dirichlet_group_t G, const acb_conrey_t a, const acb_conrey_t b);

long n_dirichlet_chi(const acb_dirichlet_group_t G, ulong m, ulong n);

void acb_dirichlet_chi_conrey(acb_t res, const acb_dirichlet_group_t G, const acb_conrey_t a, const acb_conrey_t b, slong prec);

void acb_dirichlet_chi(acb_t res, const acb_dirichlet_group_t G, ulong m, ulong n, slong prec);

/* introducing character type */

/* character = reduced exponents, keep order and number */
typedef struct
{
  ulong q;           /* modulus */
  ulong n;           /* number */
  ulong order;       /* order */
  ulong * expo;      /* reduced exponents ( order * log[k] / gcd( ) ) */
}
acb_dirichlet_char_struct;

typedef acb_dirichlet_char_struct acb_dirichlet_char_t[1];

void acb_dirichlet_char_init(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G);

void acb_dirichlet_char_clear(acb_dirichlet_char_t chi);

void acb_dirichlet_char(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G, ulong n);

long n_dirichlet_char_eval(const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, ulong n);

void acb_dirichlet_char_eval(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, ulong n, slong prec);

void n_dirichlet_char_vec(long *v, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, ulong nv);

void acb_dirichlet_char_vec(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, ulong n, slong prec);

#ifdef __cplusplus
}
#endif

#endif

