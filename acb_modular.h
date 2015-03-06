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

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#ifndef ACB_MODULAR_H
#define ACB_MODULAR_H

#include "acb.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    fmpz a;
    fmpz b;
    fmpz c;
    fmpz d;
}
psl2z_struct;

typedef psl2z_struct psl2z_t[1];

static __inline__ void
psl2z_init(psl2z_t g)
{
    fmpz_init(&g->a);
    fmpz_init(&g->b);
    fmpz_init(&g->c);
    fmpz_init(&g->d);
    fmpz_one(&g->a);
    fmpz_one(&g->d);
}

static __inline__ void
psl2z_clear(psl2z_t g)
{
    fmpz_clear(&g->a);
    fmpz_clear(&g->b);
    fmpz_clear(&g->c);
    fmpz_clear(&g->d);
}

static __inline__ void
psl2z_swap(psl2z_t f, psl2z_t g)
{
    psl2z_struct h = *f;
    *f = *g;
    *g = h;
}

static __inline__ void
psl2z_set(psl2z_t h, const psl2z_t g)
{
    fmpz_set(&h->a, &g->a);
    fmpz_set(&h->b, &g->b);
    fmpz_set(&h->c, &g->c);
    fmpz_set(&h->d, &g->d);
}

static __inline__ void
psl2z_one(psl2z_t g)
{
    fmpz_one(&g->a);
    fmpz_zero(&g->b);
    fmpz_zero(&g->c);
    fmpz_one(&g->d);
}

static __inline__ void
psl2z_print(const psl2z_t g)
{
    printf("[");
    fmpz_print(&g->a); printf(" ");
    fmpz_print(&g->b); printf("; ");
    fmpz_print(&g->c); printf(" ");
    fmpz_print(&g->d); printf("]");
}

static __inline__ int
psl2z_equal(const psl2z_t f, const psl2z_t g)
{
    return fmpz_equal(&f->a, &g->a)
        && fmpz_equal(&f->b, &g->b)
        && fmpz_equal(&f->c, &g->c)
        && fmpz_equal(&f->d, &g->d);
}

void psl2z_mul(psl2z_t h, const psl2z_t f, const psl2z_t g);

void psl2z_inv(psl2z_t h, const psl2z_t g);

int psl2z_is_one(const psl2z_t g);

int psl2z_is_correct(const psl2z_t g);

void psl2z_randtest(psl2z_t g, flint_rand_t state, long bits);

void acb_modular_transform(acb_t w, const psl2z_t g, const acb_t z, long prec);

void acb_modular_fundamental_domain_approx_d(psl2z_t g,
    double x, double y, double one_minus_eps);

void acb_modular_fundamental_domain_approx_arf(psl2z_t g,
    const arf_t xx, const arf_t yy, const arf_t one_minus_eps, long prec);

void acb_modular_fundamental_domain_approx(acb_t w, psl2z_t g, const acb_t z,
        const arf_t one_minus_eps, long prec);

int acb_modular_is_in_fundamental_domain(const acb_t z, const arf_t tol, long prec);

void acb_modular_addseq_theta(long * exponents, long * aindex, long * bindex, long num);

void acb_modular_addseq_eta(long * exponents, long * aindex, long * bindex, long num);

void acb_modular_theta_transform(int * R, int * S, int * C, const psl2z_t g);

void acb_modular_theta_sum(acb_ptr theta1, acb_ptr theta2,
        acb_ptr theta3, acb_ptr theta4,
    const acb_t w, int w_is_unit, const acb_t q, long len, long prec);

void acb_modular_theta_notransform(acb_t theta1, acb_t theta2,
    acb_t theta3, acb_t theta4, const acb_t z, const acb_t tau,
    long prec);

void acb_modular_theta(acb_t theta1, acb_t theta2,
    acb_t theta3, acb_t theta4, const acb_t z, const acb_t tau,
    long prec);

void acb_modular_j(acb_t z, const acb_t tau, long prec);

int acb_modular_epsilon_arg(const psl2z_t g);

void acb_modular_eta_sum(acb_t eta, const acb_t q, long prec);

void acb_modular_eta(acb_t z, const acb_t tau, long prec);

void acb_modular_lambda(acb_t r, const acb_t tau, long prec);

void acb_modular_delta(acb_t r, const acb_t tau, long prec);

void acb_modular_eisenstein(acb_ptr r, const acb_t tau, long len, long prec);

void acb_modular_elliptic_p(acb_t r, const acb_t z, const acb_t tau, long prec);

void acb_modular_elliptic_p_zpx(acb_ptr r, const acb_t z, const acb_t tau, long len, long prec);

void acb_modular_elliptic_k(acb_t k, const acb_t m, long prec);

void acb_modular_elliptic_k_cpx(acb_ptr w, const acb_t m, long len, long prec);

void acb_modular_elliptic_e(acb_t res, const acb_t m, long prec);

#ifdef __cplusplus
}
#endif

#endif

