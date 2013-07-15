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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#ifndef FMPRB_POLY_H
#define FMPRB_POLY_H

#include "fmprb.h"
#include "fmpcb.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

typedef struct
{
    fmprb_struct * coeffs;
    long length;
    long alloc;
}
fmprb_poly_struct;

typedef fmprb_poly_struct fmprb_poly_t[1];


/* Memory management */

void fmprb_poly_init(fmprb_poly_t poly);

void fmprb_poly_init2(fmprb_poly_t poly, long len);

void fmprb_poly_clear(fmprb_poly_t poly);

void fmprb_poly_fit_length(fmprb_poly_t poly, long len);

void _fmprb_poly_set_length(fmprb_poly_t poly, long len);

void _fmprb_poly_normalise(fmprb_poly_t poly);

static __inline__ void
fmprb_poly_swap(fmprb_poly_t poly1, fmprb_poly_t poly2)
{
    fmprb_poly_struct t = *poly1;
    *poly1 = *poly2;
    *poly2 = t;
}

void fmprb_poly_set(fmprb_poly_t poly, const fmprb_poly_t src);

/* Basic manipulation */

static __inline__ long fmprb_poly_length(const fmprb_poly_t poly)
{
    return poly->length;
}

static __inline__ long fmprb_poly_degree(const fmprb_poly_t poly)
{
    return poly->length - 1;
}

static __inline__ void fmprb_poly_zero(fmprb_poly_t poly)
{
    poly->length = 0;
}

static __inline__ void
fmprb_poly_one(fmprb_poly_t poly)
{
    fmprb_poly_fit_length(poly, 1);
    fmprb_one(poly->coeffs);
    _fmprb_poly_set_length(poly, 1);
}

void fmprb_poly_set_coeff_si(fmprb_poly_t poly, long n, long x);

void fmprb_poly_set_coeff_fmprb(fmprb_poly_t poly, long n, const fmprb_t x);

void fmprb_poly_get_coeff_fmprb(fmprb_t x, const fmprb_poly_t poly, long n);

#define fmprb_poly_get_coeff_ptr(poly, n) \
    ((n) < (poly)->length ? (poly)->coeffs + (n) : NULL)

void _fmprb_poly_reverse(fmprb_struct * res, const fmprb_struct * poly, long len, long n);

void _fmprb_poly_shift_right(fmprb_struct * res, const fmprb_struct * poly, long len, long n);

void fmprb_poly_shift_right(fmprb_poly_t res, const fmprb_poly_t poly, long n);

void _fmprb_poly_shift_left(fmprb_struct * res, const fmprb_struct * poly, long len, long n);

void fmprb_poly_shift_left(fmprb_poly_t res, const fmprb_poly_t poly, long n);

static __inline__ void
fmprb_poly_truncate(fmprb_poly_t poly, long newlen)
{
    if (poly->length > newlen)
    {
        long i;
        for (i = newlen; i < poly->length; i++)
            fmprb_zero(poly->coeffs + i);
        poly->length = newlen;
        _fmprb_poly_normalise(poly);
    }
}

/* Conversions */

void fmprb_poly_set_fmpz_poly(fmprb_poly_t poly, const fmpz_poly_t src, long prec);

void fmprb_poly_set_fmpq_poly(fmprb_poly_t poly, const fmpq_poly_t src, long prec);

static __inline__ void
fmprb_poly_set_fmprb(fmprb_poly_t poly, const fmprb_t c)
{
    fmprb_poly_fit_length(poly, 1);
    fmprb_set(poly->coeffs, c);
    _fmprb_poly_set_length(poly, !fmprb_is_zero(poly->coeffs));
}

void fmprb_poly_set_si(fmprb_poly_t poly, long c);

/* Comparisons */

int fmprb_poly_contains_fmpq_poly(const fmprb_poly_t poly1, const fmpq_poly_t poly2);

int fmprb_poly_equal(const fmprb_poly_t A, const fmprb_poly_t B);

int _fmprb_poly_overlaps(const fmprb_struct * poly1, long len1, const fmprb_struct * poly2, long len2);

int fmprb_poly_overlaps(const fmprb_poly_t poly1, const fmprb_poly_t poly2);

/* IO */

void fmprb_poly_printd(const fmprb_poly_t poly, long digits);

/* Random generation */

void fmprb_poly_randtest(fmprb_poly_t poly, flint_rand_t state, long len, long prec, long mag_bits);

/* Arithmetic */

void
_fmprb_poly_add(fmprb_struct * res, const fmprb_struct * poly1, long len1,
    const fmprb_struct * poly2, long len2, long prec);

void fmprb_poly_add(fmprb_poly_t res, const fmprb_poly_t poly1,
              const fmprb_poly_t poly2, long prec);

void _fmprb_poly_sub(fmprb_struct * res, const fmprb_struct * poly1, long len1,
    const fmprb_struct * poly2, long len2, long prec);

void fmprb_poly_sub(fmprb_poly_t res, const fmprb_poly_t poly1,
              const fmprb_poly_t poly2, long prec);

void _fmprb_poly_mullow_ztrunc(fmprb_struct * res,
    const fmprb_struct * poly1, long len1,
    const fmprb_struct * poly2, long len2, long n, long prec);

void fmprb_poly_mullow_ztrunc(fmprb_poly_t res, const fmprb_poly_t poly1,
                                            const fmprb_poly_t poly2,
                                                long n, long prec);

void _fmprb_poly_mullow_classical(fmprb_struct * res,
    const fmprb_struct * poly1, long len1,
    const fmprb_struct * poly2, long len2, long n, long prec);

void fmprb_poly_mullow_classical(fmprb_poly_t res, const fmprb_poly_t poly1,
                                            const fmprb_poly_t poly2,
                                                long n, long prec);

void _fmprb_poly_mullow_block(fmprb_struct * C,
    const fmprb_struct * A, long lenA,
    const fmprb_struct * B, long lenB, long n, long prec);

void fmprb_poly_mullow_block(fmprb_poly_t res, const fmprb_poly_t poly1,
              const fmprb_poly_t poly2, long len, long prec);

void _fmprb_poly_mullow(fmprb_struct * C,
    const fmprb_struct * A, long lenA,
    const fmprb_struct * B, long lenB, long n, long prec);

void fmprb_poly_mullow(fmprb_poly_t res, const fmprb_poly_t poly1,
              const fmprb_poly_t poly2, long len, long prec);

void _fmprb_poly_mul(fmprb_struct * C,
    const fmprb_struct * A, long lenA,
    const fmprb_struct * B, long lenB, long prec);

void fmprb_poly_mul(fmprb_poly_t res, const fmprb_poly_t poly1,
              const fmprb_poly_t poly2, long prec);

static __inline__ void
_fmprb_poly_mul_monic(fmprb_struct * res, const fmprb_struct * poly1, long len1,
    const fmprb_struct * poly2, long len2, long prec)
{
    if (len1 + len2 - 2 > 0)
        _fmprb_poly_mullow(res, poly1, len1, poly2, len2, len1 + len2 - 2, prec);
    fmprb_one(res + len1 + len2 - 2);
}

void _fmprb_poly_inv_series(fmprb_struct * Qinv,
    const fmprb_struct * Q, long Qlen, long len, long prec);

void fmprb_poly_inv_series(fmprb_poly_t Qinv, const fmprb_poly_t Q, long n, long prec);

void  _fmprb_poly_div_series(fmprb_struct * Q, const fmprb_struct * A, long Alen,
    const fmprb_struct * B, long Blen, long n, long prec);

void fmprb_poly_div_series(fmprb_poly_t Q, const fmprb_poly_t A, const fmprb_poly_t B, long n, long prec);

void
_fmprb_poly_div(fmprb_struct * Q,
    const fmprb_struct * A, long lenA,
    const fmprb_struct * B, long lenB, long prec);

void _fmprb_poly_divrem(fmprb_struct * Q, fmprb_struct * R,
    const fmprb_struct * A, long lenA,
    const fmprb_struct * B, long lenB, long prec);

void _fmprb_poly_rem(fmprb_struct * R,
    const fmprb_struct * A, long lenA,
    const fmprb_struct * B, long lenB, long prec);

void fmprb_poly_divrem(fmprb_poly_t Q, fmprb_poly_t R,
                             const fmprb_poly_t A, const fmprb_poly_t B, long prec);

void _fmprb_poly_div_root(fmprb_struct * Q, fmprb_t R, const fmprb_struct * A,
    long len, const fmprb_t c, long prec);

/* Product trees */

void _fmprb_poly_product_roots(fmprb_struct * poly, const fmprb_struct * xs, long n, long prec);

void fmprb_poly_product_roots(fmprb_poly_t poly, fmprb_struct * xs, long n, long prec);

fmprb_struct ** _fmprb_poly_tree_alloc(long len);

void _fmprb_poly_tree_free(fmprb_struct ** tree, long len);

void _fmprb_poly_tree_build(fmprb_struct ** tree, const fmprb_struct * roots, long len, long prec);

/* Composition */

void _fmprb_poly_compose(fmprb_struct * res,
    const fmprb_struct * poly1, long len1,
    const fmprb_struct * poly2, long len2, long prec);

void fmprb_poly_compose(fmprb_poly_t res,
              const fmprb_poly_t poly1, const fmprb_poly_t poly2, long prec);

void _fmprb_poly_compose_horner(fmprb_struct * res,
    const fmprb_struct * poly1, long len1,
    const fmprb_struct * poly2, long len2, long prec);

void fmprb_poly_compose_horner(fmprb_poly_t res,
              const fmprb_poly_t poly1, const fmprb_poly_t poly2, long prec);

void _fmprb_poly_compose_divconquer(fmprb_struct * res,
    const fmprb_struct * poly1, long len1,
    const fmprb_struct * poly2, long len2, long prec);

void fmprb_poly_compose_divconquer(fmprb_poly_t res,
              const fmprb_poly_t poly1, const fmprb_poly_t poly2, long prec);

void _fmprb_poly_compose_series_horner(fmprb_struct * res, const fmprb_struct * poly1, long len1,
                            const fmprb_struct * poly2, long len2, long n, long prec);

void fmprb_poly_compose_series_horner(fmprb_poly_t res,
                    const fmprb_poly_t poly1,
                    const fmprb_poly_t poly2, long n, long prec);

void _fmprb_poly_compose_series_brent_kung(fmprb_struct * res, const fmprb_struct * poly1, long len1,
                            const fmprb_struct * poly2, long len2, long n, long prec);

void fmprb_poly_compose_series_brent_kung(fmprb_poly_t res,
                    const fmprb_poly_t poly1,
                    const fmprb_poly_t poly2, long n, long prec);

void _fmprb_poly_compose_series(fmprb_struct * res, const fmprb_struct * poly1, long len1,
                            const fmprb_struct * poly2, long len2, long n, long prec);

void fmprb_poly_compose_series(fmprb_poly_t res,
                    const fmprb_poly_t poly1,
                    const fmprb_poly_t poly2, long n, long prec);

/* Reversion */

void _fmprb_poly_revert_series_lagrange(fmprb_struct * Qinv, const fmprb_struct * Q, long n, long prec);
void fmprb_poly_revert_series_lagrange(fmprb_poly_t Qinv, const fmprb_poly_t Q, long n, long prec);

void _fmprb_poly_revert_series_newton(fmprb_struct * Qinv, const fmprb_struct * Q, long n, long prec);
void fmprb_poly_revert_series_newton(fmprb_poly_t Qinv, const fmprb_poly_t Q, long n, long prec);

void _fmprb_poly_revert_series_lagrange_fast(fmprb_struct * Qinv, const fmprb_struct * Q, long n, long prec);
void fmprb_poly_revert_series_lagrange_fast(fmprb_poly_t Qinv, const fmprb_poly_t Q, long n, long prec);

void _fmprb_poly_revert_series(fmprb_struct * Qinv, const fmprb_struct * Q, long n, long prec);
void fmprb_poly_revert_series(fmprb_poly_t Qinv, const fmprb_poly_t Q, long n, long prec);

/* Evaluation and interpolation */

void _fmprb_poly_evaluate_horner(fmprb_t res, const fmprb_struct * f, long len, const fmprb_t a, long prec);
void fmprb_poly_evaluate_horner(fmprb_t res, const fmprb_poly_t f, const fmprb_t a, long prec);

void _fmprb_poly_evaluate_rectangular(fmprb_t y, const fmprb_struct * poly, long len, const fmprb_t x, long prec);
void fmprb_poly_evaluate_rectangular(fmprb_t res, const fmprb_poly_t f, const fmprb_t a, long prec);

void _fmprb_poly_evaluate(fmprb_t res, const fmprb_struct * f, long len, const fmprb_t a, long prec);
void fmprb_poly_evaluate(fmprb_t res, const fmprb_poly_t f, const fmprb_t a, long prec);

void _fmprb_poly_evaluate2_horner(fmprb_t y, fmprb_t z, const fmprb_struct * f, long len, const fmprb_t x, long prec);
void fmprb_poly_evaluate2_horner(fmprb_t y, fmprb_t z, const fmprb_poly_t f, const fmprb_t x, long prec);

void _fmprb_poly_evaluate2_rectangular(fmprb_t y, fmprb_t z, const fmprb_struct * f, long len, const fmprb_t x, long prec);
void fmprb_poly_evaluate2_rectangular(fmprb_t y, fmprb_t z, const fmprb_poly_t f, const fmprb_t x, long prec);

void _fmprb_poly_evaluate2(fmprb_t y, fmprb_t z, const fmprb_struct * f, long len, const fmprb_t x, long prec);
void fmprb_poly_evaluate2(fmprb_t y, fmprb_t z, const fmprb_poly_t f, const fmprb_t x, long prec);


void _fmprb_poly_evaluate_fmpcb_horner(fmpcb_t res, const fmprb_struct * f, long len, const fmpcb_t x, long prec);
void fmprb_poly_evaluate_fmpcb_horner(fmpcb_t res, const fmprb_poly_t f, const fmpcb_t a, long prec);

void _fmprb_poly_evaluate_fmpcb_rectangular(fmpcb_t y, const fmprb_struct * poly, long len, const fmpcb_t x, long prec);
void fmprb_poly_evaluate_fmpcb_rectangular(fmpcb_t res, const fmprb_poly_t f, const fmpcb_t a, long prec);

void _fmprb_poly_evaluate_fmpcb(fmpcb_t res, const fmprb_struct * f, long len, const fmpcb_t x, long prec);
void fmprb_poly_evaluate_fmpcb(fmpcb_t res, const fmprb_poly_t f, const fmpcb_t a, long prec);

void _fmprb_poly_evaluate2_fmpcb_horner(fmpcb_t y, fmpcb_t z, const fmprb_struct * f, long len, const fmpcb_t x, long prec);
void fmprb_poly_evaluate2_fmpcb_horner(fmpcb_t y, fmpcb_t z, const fmprb_poly_t f, const fmpcb_t x, long prec);

void _fmprb_poly_evaluate2_fmpcb_rectangular(fmpcb_t y, fmpcb_t z, const fmprb_struct * f, long len, const fmpcb_t x, long prec);
void fmprb_poly_evaluate2_fmpcb_rectangular(fmpcb_t y, fmpcb_t z, const fmprb_poly_t f, const fmpcb_t x, long prec);

void _fmprb_poly_evaluate2_fmpcb(fmpcb_t y, fmpcb_t z, const fmprb_struct * f, long len, const fmpcb_t x, long prec);
void fmprb_poly_evaluate2_fmpcb(fmpcb_t y, fmpcb_t z, const fmprb_poly_t f, const fmpcb_t x, long prec);



void _fmprb_poly_evaluate_vec_iter(fmprb_struct * ys, const fmprb_struct * poly, long plen,
    const fmprb_struct * xs, long n, long prec);

void fmprb_poly_evaluate_vec_iter(fmprb_struct * ys,
        const fmprb_poly_t poly, const fmprb_struct * xs, long n, long prec);

void _fmprb_poly_evaluate_vec_fast_precomp(fmprb_struct * vs, const fmprb_struct * poly,
    long plen, fmprb_struct ** tree, long len, long prec);

void _fmprb_poly_evaluate_vec_fast(fmprb_struct * ys, const fmprb_struct * poly, long plen,
    const fmprb_struct * xs, long n, long prec);

void fmprb_poly_evaluate_vec_fast(fmprb_struct * ys,
        const fmprb_poly_t poly, const fmprb_struct * xs, long n, long prec);

void _fmprb_poly_interpolate_newton(fmprb_struct * poly, const fmprb_struct * xs,
    const fmprb_struct * ys, long n, long prec);

void fmprb_poly_interpolate_newton(fmprb_poly_t poly,
    const fmprb_struct * xs, const fmprb_struct * ys, long n, long prec);

void
_fmprb_poly_interpolate_barycentric(fmprb_struct * poly,
    const fmprb_struct * xs, const fmprb_struct * ys, long n, long prec);

void fmprb_poly_interpolate_barycentric(fmprb_poly_t poly,
    const fmprb_struct * xs, const fmprb_struct * ys, long n, long prec);

void _fmprb_poly_interpolation_weights(fmprb_struct * w,
    fmprb_struct ** tree, long len, long prec);

void _fmprb_poly_interpolate_fast_precomp(fmprb_struct * poly,
    const fmprb_struct * ys, fmprb_struct ** tree, const fmprb_struct * weights,
    long len, long prec);

void _fmprb_poly_interpolate_fast(fmprb_struct * poly,
    const fmprb_struct * xs, const fmprb_struct * ys, long len, long prec);

void fmprb_poly_interpolate_fast(fmprb_poly_t poly,
        const fmprb_struct * xs, const fmprb_struct * ys, long n, long prec);

/* Derivative and integral */

void _fmprb_poly_derivative(fmprb_struct * res, const fmprb_struct * poly, long len, long prec);

void fmprb_poly_derivative(fmprb_poly_t res, const fmprb_poly_t poly, long prec);

void _fmprb_poly_integral(fmprb_struct * res, const fmprb_struct * poly, long len, long prec);

void fmprb_poly_integral(fmprb_poly_t res, const fmprb_poly_t poly, long prec);

/* Special functions */

void _fmprb_poly_rsqrt_series(fmprb_struct * g,
    const fmprb_struct * h, long hlen, long len, long prec);

void fmprb_poly_rsqrt_series(fmprb_poly_t g, const fmprb_poly_t h, long n, long prec);

void _fmprb_poly_sqrt_series(fmprb_struct * g,
    const fmprb_struct * h, long hlen, long len, long prec);

void fmprb_poly_sqrt_series(fmprb_poly_t g, const fmprb_poly_t h, long n, long prec);

void _fmprb_poly_log_series(fmprb_struct * res, const fmprb_struct * f, long flen, long n, long prec);

void fmprb_poly_log_series(fmprb_poly_t res, const fmprb_poly_t f, long n, long prec);

void _fmprb_poly_exp_series_basecase(fmprb_struct * f,
        const fmprb_struct * h, long hlen, long n, long prec);

void fmprb_poly_exp_series_basecase(fmprb_poly_t f, const fmprb_poly_t h, long n, long prec);

void _fmprb_poly_exp_series(fmprb_struct * f, const fmprb_struct * h, long hlen, long n, long prec);

void fmprb_poly_exp_series(fmprb_poly_t f, const fmprb_poly_t h, long n, long prec);

void fmprb_poly_log_gamma_series(fmprb_poly_t z, long n, long prec);

void _fmprb_poly_rfac_series_ui(fmprb_struct * res, const fmprb_struct * f, long flen, ulong r, long trunc, long prec);

void fmprb_poly_rfac_series_ui(fmprb_poly_t res, const fmprb_poly_t f, ulong r, long trunc, long prec);

/* Root-finding */

void _fmprb_poly_newton_convergence_factor(fmpr_t convergence_factor,
    const fmprb_struct * poly, long len,
    const fmprb_t convergence_interval, long prec);

int _fmprb_poly_newton_step(fmprb_t xnew, const fmprb_struct * poly, long len,
    const fmprb_t x,
    const fmprb_t convergence_interval,
    const fmpr_t convergence_factor, long prec);

void _fmprb_poly_newton_refine_root(fmprb_t r, const fmprb_struct * poly,
    long len, const fmprb_t start,
    const fmprb_t convergence_interval,
    const fmpr_t convergence_factor,
    long eval_extra_prec,
    long prec);

/* Macros */

#define NEWTON_ITER_BEGIN(from, to, step_from, step_to) \
    { \
        long __steps[FLINT_BITS], __i; \
        __steps[__i = 0] = step_to = (to); \
        while (step_to > (from)) \
            __steps[++__i] = (step_to = (step_to + 1) / 2); \
        for (__i--; __i >= 0; __i--) \
        { \
            step_from = step_to; \
            step_to = __steps[__i]; \


#define NEWTON_ITER_END \
        } \
    } \

#endif

