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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#ifndef ARB_CALC_H
#define ARB_CALC_H

#include "arb.h"
#include "arb_poly.h"
#include "arb_mat.h"

#ifdef __cplusplus
extern "C" {
#endif

extern TLS_PREFIX int arb_calc_verbose;

typedef int (*arb_calc_func_t)(arb_ptr out,
    const arb_t inp, void * param, long order, long prec);

#define ARB_CALC_SUCCESS 0
#define ARB_CALC_IMPRECISE_INPUT 1
#define ARB_CALC_NO_CONVERGENCE 2

/* Root-finding */

typedef struct
{
    arf_struct a;
    arf_struct b;
}
arf_interval_struct;

typedef arf_interval_struct arf_interval_t[1];
typedef arf_interval_struct * arf_interval_ptr;
typedef const arf_interval_struct * arf_interval_srcptr;

static __inline__ void
arf_interval_init(arf_interval_t v)
{
    arf_init(&v->a);
    arf_init(&v->b);
}

static __inline__ void
arf_interval_clear(arf_interval_t v)
{
    arf_clear(&v->a);
    arf_clear(&v->b);
}

static __inline__ arf_interval_ptr
_arf_interval_vec_init(long n)
{
    long i;
    arf_interval_ptr v = (arf_interval_ptr) flint_malloc(sizeof(arf_interval_struct) * n);

    for (i = 0; i < n; i++)
        arf_interval_init(v + i);

    return v;
}

static __inline__ void
_arf_interval_vec_clear(arf_interval_ptr v, long n)
{
    long i;
    for (i = 0; i < n; i++)
        arf_interval_clear(v + i);
    flint_free(v);
}

static __inline__ void
arf_interval_set(arf_interval_t v, const arf_interval_t u)
{
    arf_set(&v->a, &u->a);
    arf_set(&v->b, &u->b);
}

static __inline__ void
arf_interval_swap(arf_interval_t v, arf_interval_t u)
{
    arf_swap(&v->a, &u->a);
    arf_swap(&v->b, &u->b);
}

static __inline__ void
arf_interval_get_arb(arb_t x, const arf_interval_t v, long prec)
{
    arb_set_interval_arf(x, &v->a, &v->b, prec);
}

static __inline__ void
arf_interval_printd(const arf_interval_t v, long n)
{
    printf("[");
    arf_printd(&v->a, n);
    printf(", ");
    arf_printd(&v->b, n);
    printf("]");
}

/* bisection */

int arb_calc_partition(arf_interval_t L, arf_interval_t R,
    arb_calc_func_t func, void * param, const arf_interval_t block, long prec);

long arb_calc_isolate_roots(arf_interval_ptr * blocks, int ** flags,
    arb_calc_func_t func, void * param,
    const arf_interval_t block, long maxdepth, long maxeval, long maxfound,
    long prec);

int arb_calc_refine_root_bisect(arf_interval_t r, arb_calc_func_t func,
    void * param, const arf_interval_t start, long iter, long prec);

/* newton iteration */

void arb_calc_newton_conv_factor(arf_t conv_factor,
    arb_calc_func_t func, void * param, const arb_t conv_region, long prec);

int arb_calc_newton_step(arb_t xnew, arb_calc_func_t func, void * param,
    const arb_t x, const arb_t conv_region, const arf_t conv_factor, long prec);

int arb_calc_refine_root_newton(arb_t r, arb_calc_func_t func,
    void * param, const arb_t start, const arb_t conv_region,
    const arf_t conv_factor, long eval_extra_prec, long prec);



#ifdef __cplusplus
}
#endif

#endif

