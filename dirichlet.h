/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef DIRICHLET_H
#define DIRICHLET_H

#ifdef DIRICHLET_INLINES_C
#define DIRICHLET_INLINE
#else
#define DIRICHLET_INLINE static __inline__
#endif

#include "dlog.h"

#ifdef __cplusplus
extern "C" {
#endif

/* should this dlog pointer be in the prime or the global group? */
typedef struct
{
    ulong p;    /* underlying prime */
    int e;      /* exponent */
    nmod_t pe;  /* modulus */
    ulong phi;  /* phi(p^e) */
    ulong g;    /* conrey generator */
    dlog_precomp_struct * dlog;  /* precomputed data for discrete log mod p^e */
}
dirichlet_prime_group_struct;

typedef struct
{
    ulong q;                /* modulus */
    ulong q_even;           /* even part of modulus */
    nmod_t mod;             /* modulus with precomputed inverse */
    ulong rad_q;            /* radical = product of odd primes */
    ulong phi_q;            /* phi(q) = group size */
    slong neven;            /* number of even components (in 0,1,2)*/
    slong num;              /* number of prime components (even + odd) */
    ulong expo;             /* exponent = largest order in G */
    dirichlet_prime_group_struct * P;
    ulong * generators;     /* generators lifted mod q */
    ulong * PHI;            /* PHI(k) = expo / phi(k)        */
}
dirichlet_group_struct;

typedef dirichlet_group_struct dirichlet_group_t[1];

DIRICHLET_INLINE ulong
dirichlet_group_size(const dirichlet_group_t G)
{
    return G->phi_q;
}

void dirichlet_group_init(dirichlet_group_t G, ulong q);
void dirichlet_subgroup_init(dirichlet_group_t H, const dirichlet_group_t G, ulong h);
void dirichlet_group_clear(dirichlet_group_t G);
void dirichlet_group_dlog_precompute(dirichlet_group_t G, ulong num);
void dirichlet_group_dlog_clear(dirichlet_group_t G);

/* properties of elements without log */

ulong dirichlet_number_primitive(const dirichlet_group_t G);
ulong dirichlet_ui_conductor(const dirichlet_group_t G, ulong a);
int dirichlet_ui_parity(const dirichlet_group_t G, ulong a);
ulong dirichlet_ui_order(const dirichlet_group_t G, ulong a);

/* characters, keep both number and log */
typedef struct
{
    ulong n;           /* number */
    ulong * log;       /* s.t. prod generators[k]^log[k] = number */
}
dirichlet_char_struct;

typedef dirichlet_char_struct dirichlet_char_t[1];

void dirichlet_char_init(dirichlet_char_t x, const dirichlet_group_t G);
void dirichlet_char_clear(dirichlet_char_t x);
void dirichlet_char_print(const dirichlet_group_t G, const dirichlet_char_t x);

DIRICHLET_INLINE void
dirichlet_char_set(dirichlet_char_t x, const dirichlet_group_t G, const dirichlet_char_t y)
{
    slong k;
    x->n = y->n;
    for (k = 0; k < G->num; k++)
        x->log[k] = y->log[k];
}

DIRICHLET_INLINE int
dirichlet_char_eq(const dirichlet_char_t x, const dirichlet_char_t y)
{
    return (x->n == y->n);
}

int dirichlet_char_eq_deep(const dirichlet_group_t G, const dirichlet_char_t x, const dirichlet_char_t y);
int dirichlet_char_parity(const dirichlet_group_t G, const dirichlet_char_t x);
ulong dirichlet_char_conductor(const dirichlet_group_t G, const dirichlet_char_t x);
ulong dirichlet_char_order(const dirichlet_group_t G, const dirichlet_char_t x);

void dirichlet_char_log(dirichlet_char_t x, const dirichlet_group_t G, ulong m);
ulong dirichlet_char_exp(dirichlet_char_t x, const dirichlet_group_t G);

void dirichlet_char_index(dirichlet_char_t x, const dirichlet_group_t G, ulong j);
ulong dirichlet_index_char(const dirichlet_group_t G, const dirichlet_char_t x);

void dirichlet_char_one(dirichlet_char_t x, const dirichlet_group_t G);
void dirichlet_char_first_primitive(dirichlet_char_t x, const dirichlet_group_t G);

int dirichlet_char_next(dirichlet_char_t x, const dirichlet_group_t G);
int dirichlet_char_next_primitive(dirichlet_char_t x, const dirichlet_group_t G);

void dirichlet_char_mul(dirichlet_char_t c, const dirichlet_group_t G, const dirichlet_char_t a, const dirichlet_char_t b);
void dirichlet_char_pow(dirichlet_char_t c, const dirichlet_group_t G, const dirichlet_char_t a, ulong n);
void dirichlet_char_primitive(dirichlet_char_t y, const dirichlet_group_t G, const dirichlet_char_t x, ulong cond);

#define DIRICHLET_CHI_NULL UWORD_MAX

ulong dirichlet_ui_pairing_char(const dirichlet_group_t G, const dirichlet_char_t a, const dirichlet_char_t b);
ulong dirichlet_ui_pairing(const dirichlet_group_t G, ulong m, ulong n);

void dirichlet_pairing_char(acb_t res, const dirichlet_group_t G, const dirichlet_char_t a, const dirichlet_char_t b, slong prec);
void dirichlet_pairing(acb_t res, const dirichlet_group_t G, ulong m, ulong n, slong prec);

/* introducing character type */

/* character = reduced exponents, keep order, number and conductor */
typedef struct
{
    ulong q;           /* modulus */
    nmod_t order;       /* order */
    dirichlet_char_t x;
    ulong * expo;      /* reduced exponents ( log[k] * PHI[k] / gcd( ) ) */
    int parity;        /* 0 for even char, 1 for odd */
    ulong conductor;
}
dirichlet_fullchar_struct;

typedef dirichlet_fullchar_struct dirichlet_fullchar_t[1];

DIRICHLET_INLINE ulong
dirichlet_fullchar_order(const dirichlet_fullchar_t chi)
{
    return chi->order.n;
}

DIRICHLET_INLINE ulong
dirichlet_fullchar_conductor(const dirichlet_fullchar_t chi)
{
    return chi->conductor;
}

DIRICHLET_INLINE int
dirichlet_fullchar_parity(const dirichlet_fullchar_t chi)
{
    return chi->parity;
}

void dirichlet_fullchar_init(dirichlet_fullchar_t chi, const dirichlet_group_t G);
void dirichlet_fullchar_clear(dirichlet_fullchar_t chi);
void dirichlet_fullchar_print(const dirichlet_group_t G, const dirichlet_fullchar_t chi);

DIRICHLET_INLINE void
dirichlet_fullchar_set(dirichlet_fullchar_t chi1, const dirichlet_group_t G, const dirichlet_fullchar_t chi2)
{
    slong k;

    chi1->q = chi2->q;
    chi1->conductor = chi2->conductor;
    chi1->order = chi2->order;
    chi1->parity = chi2->parity;
    dirichlet_char_set(chi1->x, G, chi2->x);
    for (k = 0; k < G->num; k++)
        chi1->expo[k] = chi2->expo[k];
}

DIRICHLET_INLINE int
dirichlet_fullchar_eq(const dirichlet_fullchar_t chi1, const dirichlet_fullchar_t chi2)
{
    return (chi1->q == chi2->q && chi1->x->n == chi2->x->n);
}

int dirichlet_fullchar_eq_deep(const dirichlet_group_t G, const dirichlet_fullchar_t chi1, const dirichlet_fullchar_t chi2);
DIRICHLET_INLINE int
dirichlet_fullchar_is_principal(const dirichlet_fullchar_t chi)
{
    return (chi->x->n == 1);
}
DIRICHLET_INLINE int
dirichlet_fullchar_is_real(const dirichlet_fullchar_t chi)
{
    return (chi->order.n <= 2);
}

void dirichlet_fullchar(dirichlet_fullchar_t chi, const dirichlet_group_t G, ulong n);
void dirichlet_fullchar_char(dirichlet_fullchar_t chi, const dirichlet_group_t G, const dirichlet_char_t x);
void dirichlet_fullchar_set_expo(dirichlet_fullchar_t chi, const dirichlet_group_t G);
void dirichlet_fullchar_normalize(dirichlet_fullchar_t chi, const dirichlet_group_t G);
void dirichlet_fullchar_denormalize(dirichlet_fullchar_t chi, const dirichlet_group_t G);

void dirichlet_fullchar_mul(dirichlet_fullchar_t chi12, const dirichlet_group_t G, const dirichlet_fullchar_t chi1, const dirichlet_fullchar_t chi2);
void dirichlet_fullchar_primitive(dirichlet_fullchar_t chi0, const dirichlet_group_t G0, const dirichlet_group_t G, const dirichlet_fullchar_t chi);

void dirichlet_fullchar_one(dirichlet_fullchar_t chi, const dirichlet_group_t G);
void dirichlet_fullchar_first_primitive(dirichlet_fullchar_t chi, const dirichlet_group_t G);

int dirichlet_fullchar_next(dirichlet_fullchar_t chi, const dirichlet_group_t G);
int dirichlet_fullchar_next_primitive(dirichlet_fullchar_t chi, const dirichlet_group_t G);

ulong dirichlet_ui_chi_char(const dirichlet_group_t G, const dirichlet_fullchar_t chi, const dirichlet_char_t x);
ulong dirichlet_ui_chi(const dirichlet_group_t G, const dirichlet_fullchar_t chi, ulong n);

void dirichlet_ui_vec_set_null(ulong *v, const dirichlet_group_t G, slong nv);
void dirichlet_ui_chi_vec_loop(ulong *v, const dirichlet_group_t G, const dirichlet_fullchar_t chi, slong nv);
void dirichlet_ui_chi_vec_primeloop(ulong *v, const dirichlet_group_t G, const dirichlet_fullchar_t chi, slong nv);
void dirichlet_ui_chi_vec(ulong *v, const dirichlet_group_t G, const dirichlet_fullchar_t chi, slong nv);

#ifdef __cplusplus
}
#endif

#endif
