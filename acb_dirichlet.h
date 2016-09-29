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

#ifndef ACB_DIRICHLET_H
#define ACB_DIRICHLET_H

#ifdef ACB_DIRICHLET_INLINES_C
#define ACB_DIRICHLET_INLINE
#else
#define ACB_DIRICHLET_INLINE static __inline__
#endif

#include "acb.h"
#include "dlog.h"

#ifdef __cplusplus
extern "C" {
#endif

void _acb_dirichlet_euler_product_real_ui(arb_t res, ulong s,
    const signed char * chi, int mod, int reciprocal, slong prec);

void acb_dirichlet_eta(acb_t res, const acb_t s, slong prec);

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
acb_dirichlet_prime_group_struct;

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
    acb_dirichlet_prime_group_struct * P;
    ulong * generators;     /* generators lifted mod q */
    ulong * PHI;            /* PHI(k) = expo / phi(k)        */
}
acb_dirichlet_group_struct;

typedef acb_dirichlet_group_struct acb_dirichlet_group_t[1];

ACB_DIRICHLET_INLINE ulong
acb_dirichlet_group_size(const acb_dirichlet_group_t G)
{
    return G->phi_q;
}

void acb_dirichlet_group_init(acb_dirichlet_group_t G, ulong q);
void acb_dirichlet_subgroup_init(acb_dirichlet_group_t H, const acb_dirichlet_group_t G, ulong h);
void acb_dirichlet_group_clear(acb_dirichlet_group_t G);
void acb_dirichlet_group_dlog_precompute(acb_dirichlet_group_t G, ulong num);
void acb_dirichlet_group_dlog_clear(acb_dirichlet_group_t G);

/* properties of elements without log */

ulong acb_dirichlet_number_primitive(const acb_dirichlet_group_t G);
ulong acb_dirichlet_ui_conductor(const acb_dirichlet_group_t G, ulong a);
int acb_dirichlet_ui_parity(const acb_dirichlet_group_t G, ulong a);
ulong acb_dirichlet_ui_order(const acb_dirichlet_group_t G, ulong a);

/* elements of the group, keep both number and log */
typedef struct
{
    ulong n;           /* number */
    ulong * log;       /* s.t. prod generators[k]^log[k] = number */
}
acb_dirichlet_conrey_struct;

typedef acb_dirichlet_conrey_struct acb_dirichlet_conrey_t[1];

void acb_dirichlet_conrey_init(acb_dirichlet_conrey_t x, const acb_dirichlet_group_t G);
void acb_dirichlet_conrey_clear(acb_dirichlet_conrey_t x);
void acb_dirichlet_conrey_print(const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t x);

ACB_DIRICHLET_INLINE void
acb_dirichlet_conrey_set(acb_dirichlet_conrey_t x, const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t y)
{
    slong k;
    x->n = y->n;
    for (k = 0; k < G->num; k++)
        x->log[k] = y->log[k];
}

ACB_DIRICHLET_INLINE int
acb_dirichlet_conrey_eq(const acb_dirichlet_conrey_t x, const acb_dirichlet_conrey_t y)
{
    return (x->n == y->n);
}

int acb_dirichlet_conrey_eq_deep(const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t x, const acb_dirichlet_conrey_t y);
int acb_dirichlet_conrey_parity(const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t x);
ulong acb_dirichlet_conrey_conductor(const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t x);
ulong acb_dirichlet_conrey_order(const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t x);

void acb_dirichlet_conrey_log(acb_dirichlet_conrey_t x, const acb_dirichlet_group_t G, ulong m);
ulong acb_dirichlet_conrey_exp(acb_dirichlet_conrey_t x, const acb_dirichlet_group_t G);

void acb_dirichlet_conrey_index(acb_dirichlet_conrey_t x, const acb_dirichlet_group_t G, ulong j);
ulong acb_dirichlet_index_conrey(const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t x);

void acb_dirichlet_conrey_one(acb_dirichlet_conrey_t x, const acb_dirichlet_group_t G);
void acb_dirichlet_conrey_first_primitive(acb_dirichlet_conrey_t x, const acb_dirichlet_group_t G);

int acb_dirichlet_conrey_next(acb_dirichlet_conrey_t x, const acb_dirichlet_group_t G);
int acb_dirichlet_conrey_next_primitive(acb_dirichlet_conrey_t x, const acb_dirichlet_group_t G);

void acb_dirichlet_conrey_mul(acb_dirichlet_conrey_t c, const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t a, const acb_dirichlet_conrey_t b);
void acb_dirichlet_conrey_pow(acb_dirichlet_conrey_t c, const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t a, ulong n);
void acb_dirichlet_conrey_primitive(acb_dirichlet_conrey_t y, const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t x, ulong cond);

#define ACB_DIRICHLET_CHI_NULL UWORD_MAX

ulong acb_dirichlet_ui_pairing_conrey(const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t a, const acb_dirichlet_conrey_t b);
ulong acb_dirichlet_ui_pairing(const acb_dirichlet_group_t G, ulong m, ulong n);

void acb_dirichlet_pairing_conrey(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t a, const acb_dirichlet_conrey_t b, slong prec);
void acb_dirichlet_pairing(acb_t res, const acb_dirichlet_group_t G, ulong m, ulong n, slong prec);

/* introducing character type */

/* character = reduced exponents, keep order, number and conductor */
typedef struct
{
    ulong q;           /* modulus */
    nmod_t order;       /* order */
    acb_dirichlet_conrey_t x;
    ulong * expo;      /* reduced exponents ( log[k] * PHI[k] / gcd( ) ) */
    int parity;        /* 0 for even char, 1 for odd */
    ulong conductor;
}
acb_dirichlet_char_struct;

typedef acb_dirichlet_char_struct acb_dirichlet_char_t[1];

ACB_DIRICHLET_INLINE ulong
acb_dirichlet_char_order(const acb_dirichlet_char_t chi)
{
    return chi->order.n;
}

ACB_DIRICHLET_INLINE ulong
acb_dirichlet_char_conductor(const acb_dirichlet_char_t chi)
{
    return chi->conductor;
}

ACB_DIRICHLET_INLINE int
acb_dirichlet_char_parity(const acb_dirichlet_char_t chi)
{
    return chi->parity;
}

void acb_dirichlet_char_init(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G);
void acb_dirichlet_char_clear(acb_dirichlet_char_t chi);
void acb_dirichlet_char_print(const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi);

ACB_DIRICHLET_INLINE void
acb_dirichlet_char_set(acb_dirichlet_char_t chi1, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi2)
{
    slong k;

    chi1->q = chi2->q;
    chi1->conductor = chi2->conductor;
    chi1->order = chi2->order;
    chi1->parity = chi2->parity;
    acb_dirichlet_conrey_set(chi1->x, G, chi2->x);
    for (k = 0; k < G->num; k++)
        chi1->expo[k] = chi2->expo[k];
}

ACB_DIRICHLET_INLINE int
acb_dirichlet_char_eq(const acb_dirichlet_char_t chi1, const acb_dirichlet_char_t chi2)
{
    return (chi1->q == chi2->q && chi1->x->n == chi2->x->n);
}

int acb_dirichlet_char_eq_deep(const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi1, const acb_dirichlet_char_t chi2);
ACB_DIRICHLET_INLINE int
acb_dirichlet_char_is_principal(const acb_dirichlet_char_t chi)
{
    return (chi->x->n == 1);
}
ACB_DIRICHLET_INLINE int
acb_dirichlet_char_is_real(const acb_dirichlet_char_t chi)
{
    return (chi->order.n <= 2);
}

void acb_dirichlet_char(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G, ulong n);
void acb_dirichlet_char_conrey(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t x);
void acb_dirichlet_char_set_expo(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G);
void acb_dirichlet_char_normalize(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G);
void acb_dirichlet_char_denormalize(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G);

void acb_dirichlet_char_mul(acb_dirichlet_char_t chi12, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi1, const acb_dirichlet_char_t chi2);
void acb_dirichlet_char_primitive(acb_dirichlet_char_t chi0, const acb_dirichlet_group_t G0, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi);

void acb_dirichlet_char_one(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G);
void acb_dirichlet_char_first_primitive(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G);

int acb_dirichlet_char_next(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G);
int acb_dirichlet_char_next_primitive(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G);

ulong acb_dirichlet_ui_chi_conrey(const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, const acb_dirichlet_conrey_t x);
ulong acb_dirichlet_ui_chi(const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, ulong n);

void acb_dirichlet_ui_vec_set_null(ulong *v, const acb_dirichlet_group_t G, slong nv);
void acb_dirichlet_ui_chi_vec_loop(ulong *v, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong nv);
void acb_dirichlet_ui_chi_vec_primeloop(ulong *v, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong nv);
void acb_dirichlet_ui_chi_vec(ulong *v, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong nv);

/* precompute powers of a root of unity */
typedef struct
{
    ulong order;
    acb_t z;
    slong size;
    slong depth;
    acb_ptr * Z;
}
acb_dirichlet_powers_struct;

typedef acb_dirichlet_powers_struct acb_dirichlet_powers_t[1];

void _acb_dirichlet_powers_init(acb_dirichlet_powers_t t, ulong order, slong size, slong depth, slong prec);
void acb_dirichlet_powers_init(acb_dirichlet_powers_t t, ulong order, slong num, slong prec);
void acb_dirichlet_powers_clear(acb_dirichlet_powers_t t);
void acb_dirichlet_power(acb_t z, const acb_dirichlet_powers_t t, ulong n, slong prec);

void acb_dirichlet_nth_root(acb_t res, ulong order, slong prec);
void acb_dirichlet_vec_nth_roots(acb_ptr z, slong len, slong prec);

void acb_dirichlet_chi(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, ulong n, slong prec);
void acb_dirichlet_chi_vec(acb_ptr v, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong nv, slong prec);

void acb_dirichlet_arb_quadratic_powers(arb_ptr v, slong nv, const arb_t x, slong prec);
void acb_dirichlet_qseries_arb(acb_t res, acb_srcptr a, const arb_t x, slong len, slong prec);
void acb_dirichlet_qseries_arb_powers_naive(acb_t res, const arb_t x, int parity, const ulong *a, const acb_dirichlet_powers_t z, slong len, slong prec);
void acb_dirichlet_qseries_arb_powers_smallorder(acb_t res, const arb_t x, int parity, const ulong *a, const acb_dirichlet_powers_t z, slong len, slong prec);

ulong acb_dirichlet_theta_length_d(ulong q, double x, slong prec);
ulong acb_dirichlet_theta_length(ulong q, const arb_t x, slong prec);
void mag_tail_kexpk2_arb(mag_t res, const arb_t a, ulong n);

void _acb_dirichlet_theta_argument_at_arb(arb_t xt, ulong q, const arb_t t, slong prec);
void acb_dirichlet_theta_arb(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, const arb_t t, slong prec);
void acb_dirichlet_ui_theta_arb(acb_t res, const acb_dirichlet_group_t G, ulong a, const arb_t t, slong prec);

void acb_dirichlet_gauss_sum_naive(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong prec);
void acb_dirichlet_gauss_sum_factor(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong prec);
void acb_dirichlet_gauss_sum_order2(acb_t res, const acb_dirichlet_char_t chi, slong prec);
void acb_dirichlet_gauss_sum_theta(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong prec);
void acb_dirichlet_gauss_sum(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong prec);

void acb_dirichlet_root_number_theta(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong prec);
void acb_dirichlet_root_number(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong prec);

void acb_dirichlet_si_poly_evaluate(acb_t res, slong * v, slong len, const acb_t z, slong prec);

void acb_dirichlet_jacobi_sum_naive(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi1, const acb_dirichlet_char_t chi2, slong prec);
ulong jacobi_one_prime(ulong p, ulong e, ulong pe, ulong cond);
void acb_dirichlet_jacobi_sum_factor(acb_t res,  const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi1, const acb_dirichlet_char_t chi2, slong prec);
void acb_dirichlet_jacobi_sum_gauss(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi1, const acb_dirichlet_char_t chi2, slong prec);
void acb_dirichlet_jacobi_sum(acb_t res, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi1, const acb_dirichlet_char_t chi2, slong prec);
void acb_dirichlet_jacobi_sum_ui(acb_t res, const acb_dirichlet_group_t G, ulong a, ulong b, slong prec);

void acb_dirichlet_l_hurwitz(acb_t res, const acb_t s, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong prec);
void acb_dirichlet_l_vec_hurwitz(acb_ptr res, const acb_t s, const acb_dirichlet_group_t G, slong prec);

/* Discrete Fourier Transform */

void _acb_dirichlet_dft_pol(acb_ptr w, acb_srcptr v, slong dv, acb_srcptr z, slong dz, slong len, slong prec);
void acb_dirichlet_dft_pol(acb_ptr w, acb_srcptr v, slong len, slong prec);
void acb_dirichlet_dft_crt(acb_ptr w, acb_srcptr v, slong len, slong prec);
void acb_dirichlet_dft_cyc(acb_ptr w, acb_srcptr v, slong len, slong prec);
void acb_dirichlet_dft_prod(acb_ptr w, acb_srcptr v, slong * cyc, slong num, slong prec);

void acb_dirichlet_dft_conrey(acb_ptr w, acb_srcptr v, const acb_dirichlet_group_t G, slong prec);
void acb_dirichlet_dft(acb_ptr w, acb_srcptr v, const acb_dirichlet_group_t G, slong prec);

#define CRT_MAX 15
typedef struct
{
    int num;
    nmod_t n;
    slong m[CRT_MAX];
    ulong M[CRT_MAX];
    ulong vM[CRT_MAX];
}
crt_struct;

typedef crt_struct crt_t[1];

void crt_init(crt_t c, ulong n);

void crt_decomp(acb_ptr y, acb_srcptr x, const crt_t c, ulong len);
void crt_recomp(acb_ptr y, acb_srcptr x, const crt_t c, ulong len);

typedef struct acb_dirichlet_dft_step_struct acb_dirichlet_dft_step_struct;
typedef acb_dirichlet_dft_step_struct * acb_dirichlet_dft_step_ptr;

typedef struct
{
    slong n;
    int zclear;
    acb_ptr z;
    slong num;
    acb_dirichlet_dft_step_ptr cyc;
}
acb_dirichlet_dft_cyc_struct;

typedef acb_dirichlet_dft_cyc_struct acb_dirichlet_dft_cyc_t[1];

typedef struct
{
    slong n;
    slong num;
    acb_dirichlet_dft_step_ptr cyc;
}
acb_dirichlet_dft_prod_struct;

typedef acb_dirichlet_dft_prod_struct acb_dirichlet_dft_prod_t[1];

typedef struct
{
    slong n;
    crt_t c;
    /* then a product */
    acb_dirichlet_dft_step_ptr cyc;
}
acb_dirichlet_dft_crt_struct;

typedef acb_dirichlet_dft_crt_struct acb_dirichlet_dft_crt_t[1];

typedef struct
{
    slong n;
    slong dv;
    int zclear;
    acb_ptr z;
    slong dz;
}
acb_dirichlet_dft_pol_struct;

typedef acb_dirichlet_dft_pol_struct acb_dirichlet_dft_pol_t[1];

typedef struct
{
    int type;
    union
    {
        acb_dirichlet_dft_cyc_t cyc;
        acb_dirichlet_dft_prod_t prod;
        acb_dirichlet_dft_crt_t crt;
        acb_dirichlet_dft_pol_t pol;
    } t;
}
acb_dirichlet_dft_pre_struct;

typedef acb_dirichlet_dft_pre_struct acb_dirichlet_dft_pre_t[1];

/* covers both product and cyclic case */
struct
acb_dirichlet_dft_step_struct
{
    /* [G:H] */
    slong m;
    /* card H */
    slong M;
    slong dv; /* = M for prod, also = M if cyc is reordered */
    /* pointer on some roots of unity, if needed */
    acb_srcptr z;
    /* index of mM in z */
    slong dz;
    /* dft to call on H */
    acb_dirichlet_dft_pre_t pre;
    /* dft to call on G/H ? */
};

/*typedef acb_dirichlet_dft_pre_struct acb_dirichlet_dft_pre_t[1];*/

enum
{
    DFT_POL, DFT_CYC, DFT_PROD, DFT_CRT /*, DFT_2E, DFT_CONV */
};

void acb_dirichlet_dft_step(acb_ptr w, acb_srcptr v, acb_dirichlet_dft_step_ptr cyc, slong num, slong prec);

void acb_dirichlet_dft_precomp(acb_ptr w, acb_srcptr v, const acb_dirichlet_dft_pre_t pre, slong prec);

ACB_DIRICHLET_INLINE void
acb_dirichlet_dft_pol_precomp(acb_ptr w, acb_srcptr v, const acb_dirichlet_dft_pol_t pol, slong prec)
{
    _acb_dirichlet_dft_pol(w, v, pol->dv, pol->z, pol->dz, pol->n, prec);
}
ACB_DIRICHLET_INLINE void
acb_dirichlet_dft_cyc_precomp(acb_ptr w, acb_srcptr v, const acb_dirichlet_dft_cyc_t cyc, slong prec)
{
    acb_dirichlet_dft_step(w, v, cyc->cyc, cyc->num, prec);
}
void acb_dirichlet_dft_crt_precomp(acb_ptr w, acb_srcptr v, const acb_dirichlet_dft_crt_t crt, slong prec);
void acb_dirichlet_dft_prod_precomp(acb_ptr w, acb_srcptr v, const acb_dirichlet_dft_prod_t prod, slong prec);

void _acb_dirichlet_dft_precomp_init(acb_dirichlet_dft_pre_t pre, slong dv, acb_ptr z, slong dz, slong len, slong prec);
void acb_dirichlet_dft_precomp_init(acb_dirichlet_dft_pre_t pre, slong len, slong prec);
void acb_dirichlet_dft_precomp_clear(acb_dirichlet_dft_pre_t pre);

acb_dirichlet_dft_step_ptr _acb_dirichlet_dft_steps_prod(slong * m, slong num, slong prec);

ACB_DIRICHLET_INLINE void
acb_dirichlet_dft_prod_init(acb_dirichlet_dft_prod_t t, slong * cyc, slong num, slong prec)
{
    t->num = num;
    t->cyc = _acb_dirichlet_dft_steps_prod(cyc, num, prec);
}

void acb_dirichlet_dft_prod_clear(acb_dirichlet_dft_prod_t t);

void _acb_dirichlet_dft_cyc_init_z_fac(acb_dirichlet_dft_cyc_t t, n_factor_t fac, slong dv, acb_ptr z, slong dz, slong len, slong prec);
void _acb_dirichlet_dft_cyc_init(acb_dirichlet_dft_cyc_t t, slong dv, slong len, slong prec);

ACB_DIRICHLET_INLINE void
acb_dirichlet_dft_cyc_init(acb_dirichlet_dft_cyc_t t, slong len, slong prec)
{
    _acb_dirichlet_dft_cyc_init(t, 1, len, prec);
}

void acb_dirichlet_dft_cyc_clear(acb_dirichlet_dft_cyc_t t);

void _acb_dirichlet_dft_pol_init(acb_dirichlet_dft_pol_t pol, slong dv, acb_ptr z, slong dz, slong len, slong prec);

ACB_DIRICHLET_INLINE void
acb_dirichlet_dft_pol_init(acb_dirichlet_dft_pol_t pol, slong len, slong prec)
{
    _acb_dirichlet_dft_pol_init(pol, 1, NULL, 0, len, prec);
}

ACB_DIRICHLET_INLINE void
acb_dirichlet_dft_pol_clear(acb_dirichlet_dft_pol_t pol)
{
    if (pol->zclear)
    {
        flint_printf("  ## clearing pol [len=%ld]....", pol->n);
        _acb_vec_clear(pol->z, pol->n);
        flint_printf("done\n");
    }
}

ACB_DIRICHLET_INLINE void
acb_dirichlet_dft_crt_init(acb_dirichlet_dft_crt_t crt, slong len, slong prec)
{
    crt->n = len;
    crt_init(crt->c, len);
    crt->cyc = _acb_dirichlet_dft_steps_prod(crt->c->m, crt->c->num, prec);
}

ACB_DIRICHLET_INLINE void
acb_dirichlet_dft_crt_clear(acb_dirichlet_dft_crt_t crt)
{
    flint_free(crt->cyc);
}

/* utils */

ACB_DIRICHLET_INLINE void
acb_vec_printd(acb_srcptr vec, slong len, slong digits)
{
    slong i;
    for (i = 0; i < len; i++)
        acb_printd(vec + i, digits), flint_printf("\n");
}

ACB_DIRICHLET_INLINE void
acb_vec_printd_index(acb_srcptr vec, slong len, slong digits)
{
    slong i;
    for (i = 0; i < len; i++)
        flint_printf("[%ld] ",i), acb_printd(vec + i, digits), flint_printf("\n");
}

#ifdef __cplusplus
}
#endif

#endif
