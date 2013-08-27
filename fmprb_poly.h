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
    fmprb_ptr coeffs;
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

void _fmprb_poly_reverse(fmprb_ptr res, fmprb_srcptr poly, long len, long n);

void _fmprb_poly_shift_right(fmprb_ptr res, fmprb_srcptr poly, long len, long n);

void fmprb_poly_shift_right(fmprb_poly_t res, const fmprb_poly_t poly, long n);

void _fmprb_poly_shift_left(fmprb_ptr res, fmprb_srcptr poly, long len, long n);

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

int fmprb_poly_contains(const fmprb_poly_t poly1, const fmprb_poly_t poly2);

int fmprb_poly_contains_fmpz_poly(const fmprb_poly_t poly1, const fmpz_poly_t poly2);

int fmprb_poly_contains_fmpq_poly(const fmprb_poly_t poly1, const fmpq_poly_t poly2);

int fmprb_poly_equal(const fmprb_poly_t A, const fmprb_poly_t B);

int _fmprb_poly_overlaps(fmprb_srcptr poly1, long len1, fmprb_srcptr poly2, long len2);

int fmprb_poly_overlaps(const fmprb_poly_t poly1, const fmprb_poly_t poly2);

/* IO */

void fmprb_poly_printd(const fmprb_poly_t poly, long digits);

/* Random generation */

void fmprb_poly_randtest(fmprb_poly_t poly, flint_rand_t state, long len, long prec, long mag_bits);

/* Arithmetic */

void
_fmprb_poly_add(fmprb_ptr res, fmprb_srcptr poly1, long len1,
    fmprb_srcptr poly2, long len2, long prec);

void fmprb_poly_add(fmprb_poly_t res, const fmprb_poly_t poly1,
              const fmprb_poly_t poly2, long prec);

void _fmprb_poly_sub(fmprb_ptr res, fmprb_srcptr poly1, long len1,
    fmprb_srcptr poly2, long len2, long prec);

void fmprb_poly_sub(fmprb_poly_t res, const fmprb_poly_t poly1,
              const fmprb_poly_t poly2, long prec);

static __inline__ void
fmprb_poly_neg(fmprb_poly_t res, const fmprb_poly_t poly)
{
    fmprb_poly_fit_length(res, poly->length);
    _fmprb_vec_neg(res->coeffs, poly->coeffs, poly->length);
    _fmprb_poly_set_length(res, poly->length);
}

static __inline__ void
fmprb_poly_scalar_mul_2exp_si(fmprb_poly_t res, const fmprb_poly_t poly, long c)
{
    fmprb_poly_fit_length(res, poly->length);
    _fmprb_vec_scalar_mul_2exp_si(res->coeffs, poly->coeffs, poly->length, c);
    _fmprb_poly_set_length(res, poly->length);
}

void _fmprb_poly_mullow_ztrunc(fmprb_ptr res,
    fmprb_srcptr poly1, long len1,
    fmprb_srcptr poly2, long len2, long n, long prec);

void fmprb_poly_mullow_ztrunc(fmprb_poly_t res, const fmprb_poly_t poly1,
                                            const fmprb_poly_t poly2,
                                                long n, long prec);

void _fmprb_poly_mullow_classical(fmprb_ptr res,
    fmprb_srcptr poly1, long len1,
    fmprb_srcptr poly2, long len2, long n, long prec);

void fmprb_poly_mullow_classical(fmprb_poly_t res, const fmprb_poly_t poly1,
                                            const fmprb_poly_t poly2,
                                                long n, long prec);

void _fmprb_poly_mullow_block(fmprb_ptr C,
    fmprb_srcptr A, long lenA,
    fmprb_srcptr B, long lenB, long n, long prec);

void fmprb_poly_mullow_block(fmprb_poly_t res, const fmprb_poly_t poly1,
              const fmprb_poly_t poly2, long len, long prec);

void _fmprb_poly_mullow_block_scaled(fmprb_ptr z, fmprb_srcptr x, long xlen, fmprb_srcptr y, long ylen, long len, long prec);

void fmprb_poly_mullow_block_scaled(fmprb_poly_t res, const fmprb_poly_t poly1, const fmprb_poly_t poly2, long len, long prec);

void _fmprb_poly_mullow(fmprb_ptr C,
    fmprb_srcptr A, long lenA,
    fmprb_srcptr B, long lenB, long n, long prec);

void fmprb_poly_mullow(fmprb_poly_t res, const fmprb_poly_t poly1,
              const fmprb_poly_t poly2, long len, long prec);

void _fmprb_poly_mul(fmprb_ptr C,
    fmprb_srcptr A, long lenA,
    fmprb_srcptr B, long lenB, long prec);

void fmprb_poly_mul(fmprb_poly_t res, const fmprb_poly_t poly1,
              const fmprb_poly_t poly2, long prec);

static __inline__ void
_fmprb_poly_mul_monic(fmprb_ptr res, fmprb_srcptr poly1, long len1,
    fmprb_srcptr poly2, long len2, long prec)
{
    if (len1 + len2 - 2 > 0)
        _fmprb_poly_mullow(res, poly1, len1, poly2, len2, len1 + len2 - 2, prec);
    fmprb_one(res + len1 + len2 - 2);
}

void _fmprb_poly_inv_series(fmprb_ptr Qinv,
    fmprb_srcptr Q, long Qlen, long len, long prec);

void fmprb_poly_inv_series(fmprb_poly_t Qinv, const fmprb_poly_t Q, long n, long prec);

void  _fmprb_poly_div_series(fmprb_ptr Q, fmprb_srcptr A, long Alen,
    fmprb_srcptr B, long Blen, long n, long prec);

void fmprb_poly_div_series(fmprb_poly_t Q, const fmprb_poly_t A, const fmprb_poly_t B, long n, long prec);

void
_fmprb_poly_div(fmprb_ptr Q,
    fmprb_srcptr A, long lenA,
    fmprb_srcptr B, long lenB, long prec);

void _fmprb_poly_divrem(fmprb_ptr Q, fmprb_ptr R,
    fmprb_srcptr A, long lenA,
    fmprb_srcptr B, long lenB, long prec);

void _fmprb_poly_rem(fmprb_ptr R,
    fmprb_srcptr A, long lenA,
    fmprb_srcptr B, long lenB, long prec);

void fmprb_poly_divrem(fmprb_poly_t Q, fmprb_poly_t R,
                             const fmprb_poly_t A, const fmprb_poly_t B, long prec);

void _fmprb_poly_div_root(fmprb_ptr Q, fmprb_t R, fmprb_srcptr A,
    long len, const fmprb_t c, long prec);

/* Product trees */

void _fmprb_poly_product_roots(fmprb_ptr poly, fmprb_srcptr xs, long n, long prec);

void fmprb_poly_product_roots(fmprb_poly_t poly, fmprb_ptr xs, long n, long prec);

fmprb_ptr * _fmprb_poly_tree_alloc(long len);

void _fmprb_poly_tree_free(fmprb_ptr * tree, long len);

void _fmprb_poly_tree_build(fmprb_ptr * tree, fmprb_srcptr roots, long len, long prec);

/* Composition */

void _fmprb_poly_compose(fmprb_ptr res,
    fmprb_srcptr poly1, long len1,
    fmprb_srcptr poly2, long len2, long prec);

void fmprb_poly_compose(fmprb_poly_t res,
              const fmprb_poly_t poly1, const fmprb_poly_t poly2, long prec);

void _fmprb_poly_compose_horner(fmprb_ptr res,
    fmprb_srcptr poly1, long len1,
    fmprb_srcptr poly2, long len2, long prec);

void fmprb_poly_compose_horner(fmprb_poly_t res,
              const fmprb_poly_t poly1, const fmprb_poly_t poly2, long prec);

void _fmprb_poly_compose_divconquer(fmprb_ptr res,
    fmprb_srcptr poly1, long len1,
    fmprb_srcptr poly2, long len2, long prec);

void fmprb_poly_compose_divconquer(fmprb_poly_t res,
              const fmprb_poly_t poly1, const fmprb_poly_t poly2, long prec);

void _fmprb_poly_compose_series_horner(fmprb_ptr res, fmprb_srcptr poly1, long len1,
                            fmprb_srcptr poly2, long len2, long n, long prec);

void fmprb_poly_compose_series_horner(fmprb_poly_t res,
                    const fmprb_poly_t poly1,
                    const fmprb_poly_t poly2, long n, long prec);

void _fmprb_poly_compose_series_brent_kung(fmprb_ptr res, fmprb_srcptr poly1, long len1,
                            fmprb_srcptr poly2, long len2, long n, long prec);

void fmprb_poly_compose_series_brent_kung(fmprb_poly_t res,
                    const fmprb_poly_t poly1,
                    const fmprb_poly_t poly2, long n, long prec);

void _fmprb_poly_compose_series(fmprb_ptr res, fmprb_srcptr poly1, long len1,
                            fmprb_srcptr poly2, long len2, long n, long prec);

void fmprb_poly_compose_series(fmprb_poly_t res,
                    const fmprb_poly_t poly1,
                    const fmprb_poly_t poly2, long n, long prec);

/* Reversion */

void _fmprb_poly_revert_series_lagrange(fmprb_ptr Qinv, fmprb_srcptr Q, long n, long prec);
void fmprb_poly_revert_series_lagrange(fmprb_poly_t Qinv, const fmprb_poly_t Q, long n, long prec);

void _fmprb_poly_revert_series_newton(fmprb_ptr Qinv, fmprb_srcptr Q, long n, long prec);
void fmprb_poly_revert_series_newton(fmprb_poly_t Qinv, const fmprb_poly_t Q, long n, long prec);

void _fmprb_poly_revert_series_lagrange_fast(fmprb_ptr Qinv, fmprb_srcptr Q, long n, long prec);
void fmprb_poly_revert_series_lagrange_fast(fmprb_poly_t Qinv, const fmprb_poly_t Q, long n, long prec);

void _fmprb_poly_revert_series(fmprb_ptr Qinv, fmprb_srcptr Q, long n, long prec);
void fmprb_poly_revert_series(fmprb_poly_t Qinv, const fmprb_poly_t Q, long n, long prec);

/* Evaluation and interpolation */

void _fmprb_poly_evaluate_horner(fmprb_t res, fmprb_srcptr f, long len, const fmprb_t a, long prec);
void fmprb_poly_evaluate_horner(fmprb_t res, const fmprb_poly_t f, const fmprb_t a, long prec);

void _fmprb_poly_evaluate_rectangular(fmprb_t y, fmprb_srcptr poly, long len, const fmprb_t x, long prec);
void fmprb_poly_evaluate_rectangular(fmprb_t res, const fmprb_poly_t f, const fmprb_t a, long prec);

void _fmprb_poly_evaluate(fmprb_t res, fmprb_srcptr f, long len, const fmprb_t a, long prec);
void fmprb_poly_evaluate(fmprb_t res, const fmprb_poly_t f, const fmprb_t a, long prec);

void _fmprb_poly_evaluate2_horner(fmprb_t y, fmprb_t z, fmprb_srcptr f, long len, const fmprb_t x, long prec);
void fmprb_poly_evaluate2_horner(fmprb_t y, fmprb_t z, const fmprb_poly_t f, const fmprb_t x, long prec);

void _fmprb_poly_evaluate2_rectangular(fmprb_t y, fmprb_t z, fmprb_srcptr f, long len, const fmprb_t x, long prec);
void fmprb_poly_evaluate2_rectangular(fmprb_t y, fmprb_t z, const fmprb_poly_t f, const fmprb_t x, long prec);

void _fmprb_poly_evaluate2(fmprb_t y, fmprb_t z, fmprb_srcptr f, long len, const fmprb_t x, long prec);
void fmprb_poly_evaluate2(fmprb_t y, fmprb_t z, const fmprb_poly_t f, const fmprb_t x, long prec);


void _fmprb_poly_evaluate_fmpcb_horner(fmpcb_t res, fmprb_srcptr f, long len, const fmpcb_t x, long prec);
void fmprb_poly_evaluate_fmpcb_horner(fmpcb_t res, const fmprb_poly_t f, const fmpcb_t a, long prec);

void _fmprb_poly_evaluate_fmpcb_rectangular(fmpcb_t y, fmprb_srcptr poly, long len, const fmpcb_t x, long prec);
void fmprb_poly_evaluate_fmpcb_rectangular(fmpcb_t res, const fmprb_poly_t f, const fmpcb_t a, long prec);

void _fmprb_poly_evaluate_fmpcb(fmpcb_t res, fmprb_srcptr f, long len, const fmpcb_t x, long prec);
void fmprb_poly_evaluate_fmpcb(fmpcb_t res, const fmprb_poly_t f, const fmpcb_t a, long prec);

void _fmprb_poly_evaluate2_fmpcb_horner(fmpcb_t y, fmpcb_t z, fmprb_srcptr f, long len, const fmpcb_t x, long prec);
void fmprb_poly_evaluate2_fmpcb_horner(fmpcb_t y, fmpcb_t z, const fmprb_poly_t f, const fmpcb_t x, long prec);

void _fmprb_poly_evaluate2_fmpcb_rectangular(fmpcb_t y, fmpcb_t z, fmprb_srcptr f, long len, const fmpcb_t x, long prec);
void fmprb_poly_evaluate2_fmpcb_rectangular(fmpcb_t y, fmpcb_t z, const fmprb_poly_t f, const fmpcb_t x, long prec);

void _fmprb_poly_evaluate2_fmpcb(fmpcb_t y, fmpcb_t z, fmprb_srcptr f, long len, const fmpcb_t x, long prec);
void fmprb_poly_evaluate2_fmpcb(fmpcb_t y, fmpcb_t z, const fmprb_poly_t f, const fmpcb_t x, long prec);



void _fmprb_poly_evaluate_vec_iter(fmprb_ptr ys, fmprb_srcptr poly, long plen,
    fmprb_srcptr xs, long n, long prec);

void fmprb_poly_evaluate_vec_iter(fmprb_ptr ys,
        const fmprb_poly_t poly, fmprb_srcptr xs, long n, long prec);

void _fmprb_poly_evaluate_vec_fast_precomp(fmprb_ptr vs, fmprb_srcptr poly,
    long plen, fmprb_ptr * tree, long len, long prec);

void _fmprb_poly_evaluate_vec_fast(fmprb_ptr ys, fmprb_srcptr poly, long plen,
    fmprb_srcptr xs, long n, long prec);

void fmprb_poly_evaluate_vec_fast(fmprb_ptr ys,
        const fmprb_poly_t poly, fmprb_srcptr xs, long n, long prec);

void _fmprb_poly_interpolate_newton(fmprb_ptr poly, fmprb_srcptr xs,
    fmprb_srcptr ys, long n, long prec);

void fmprb_poly_interpolate_newton(fmprb_poly_t poly,
    fmprb_srcptr xs, fmprb_srcptr ys, long n, long prec);

void
_fmprb_poly_interpolate_barycentric(fmprb_ptr poly,
    fmprb_srcptr xs, fmprb_srcptr ys, long n, long prec);

void fmprb_poly_interpolate_barycentric(fmprb_poly_t poly,
    fmprb_srcptr xs, fmprb_srcptr ys, long n, long prec);

void _fmprb_poly_interpolation_weights(fmprb_ptr w,
    fmprb_ptr * tree, long len, long prec);

void _fmprb_poly_interpolate_fast_precomp(fmprb_ptr poly,
    fmprb_srcptr ys, fmprb_ptr * tree, fmprb_srcptr weights,
    long len, long prec);

void _fmprb_poly_interpolate_fast(fmprb_ptr poly,
    fmprb_srcptr xs, fmprb_srcptr ys, long len, long prec);

void fmprb_poly_interpolate_fast(fmprb_poly_t poly,
        fmprb_srcptr xs, fmprb_srcptr ys, long n, long prec);

/* Derivative and integral */

void _fmprb_poly_derivative(fmprb_ptr res, fmprb_srcptr poly, long len, long prec);

void fmprb_poly_derivative(fmprb_poly_t res, const fmprb_poly_t poly, long prec);

void _fmprb_poly_integral(fmprb_ptr res, fmprb_srcptr poly, long len, long prec);

void fmprb_poly_integral(fmprb_poly_t res, const fmprb_poly_t poly, long prec);

/* Transforms */

void fmprb_poly_borel_transform(fmprb_poly_t res, const fmprb_poly_t poly, long prec);

void _fmprb_poly_borel_transform(fmprb_ptr res, fmprb_srcptr poly, long len, long prec);

void fmprb_poly_inv_borel_transform(fmprb_poly_t res, const fmprb_poly_t poly, long prec);

void _fmprb_poly_inv_borel_transform(fmprb_ptr res, fmprb_srcptr poly, long len, long prec);

void _fmprb_poly_binomial_transform_basecase(fmprb_ptr b, fmprb_srcptr a, long alen, long len, long prec);

void fmprb_poly_binomial_transform_basecase(fmprb_poly_t b, const fmprb_poly_t a, long len, long prec);

void _fmprb_poly_binomial_transform_convolution(fmprb_ptr b, fmprb_srcptr a, long alen, long len, long prec);

void fmprb_poly_binomial_transform_convolution(fmprb_poly_t b, const fmprb_poly_t a, long len, long prec);

void _fmprb_poly_binomial_transform(fmprb_ptr b, fmprb_srcptr a, long alen, long len, long prec);

void fmprb_poly_binomial_transform(fmprb_poly_t b, const fmprb_poly_t a, long len, long prec);

/* Special functions */

void _fmprb_poly_pow_ui_trunc_binexp(fmprb_ptr res,
    fmprb_srcptr f, long flen, ulong exp, long len, long prec);

void fmprb_poly_pow_ui_trunc_binexp(fmprb_poly_t res,
    const fmprb_poly_t poly, ulong exp, long len, long prec);

void _fmprb_poly_pow_ui(fmprb_ptr res, fmprb_srcptr f, long flen, ulong exp, long prec);

void fmprb_poly_pow_ui(fmprb_poly_t res, const fmprb_poly_t poly, ulong exp, long prec);

void _fmprb_poly_rsqrt_series(fmprb_ptr g,
    fmprb_srcptr h, long hlen, long len, long prec);

void fmprb_poly_rsqrt_series(fmprb_poly_t g, const fmprb_poly_t h, long n, long prec);

void _fmprb_poly_sqrt_series(fmprb_ptr g,
    fmprb_srcptr h, long hlen, long len, long prec);

void fmprb_poly_sqrt_series(fmprb_poly_t g, const fmprb_poly_t h, long n, long prec);

void _fmprb_poly_log_series(fmprb_ptr res, fmprb_srcptr f, long flen, long n, long prec);

void fmprb_poly_log_series(fmprb_poly_t res, const fmprb_poly_t f, long n, long prec);

void _fmprb_poly_atan_series(fmprb_ptr res, fmprb_srcptr f, long flen, long n, long prec);

void fmprb_poly_atan_series(fmprb_poly_t res, const fmprb_poly_t f, long n, long prec);

void _fmprb_poly_asin_series(fmprb_ptr res, fmprb_srcptr f, long flen, long n, long prec);

void fmprb_poly_asin_series(fmprb_poly_t res, const fmprb_poly_t f, long n, long prec);

void _fmprb_poly_acos_series(fmprb_ptr res, fmprb_srcptr f, long flen, long n, long prec);

void fmprb_poly_acos_series(fmprb_poly_t res, const fmprb_poly_t f, long n, long prec);

void _fmprb_poly_exp_series_basecase(fmprb_ptr f,
        fmprb_srcptr h, long hlen, long n, long prec);

void fmprb_poly_exp_series_basecase(fmprb_poly_t f, const fmprb_poly_t h, long n, long prec);

void _fmprb_poly_exp_series(fmprb_ptr f, fmprb_srcptr h, long hlen, long n, long prec);

void fmprb_poly_exp_series(fmprb_poly_t f, const fmprb_poly_t h, long n, long prec);

void _fmprb_poly_sin_cos_series_basecase(fmprb_ptr s,
                                    fmprb_ptr c, fmprb_srcptr h, long hlen, long n, long prec);

void fmprb_poly_sin_cos_series_basecase(fmprb_poly_t s, fmprb_poly_t c,
        const fmprb_poly_t h, long n, long prec);

void _fmprb_poly_sin_cos_series_tangent(fmprb_ptr s, fmprb_ptr c,
                        const fmprb_srcptr h, long hlen, long len, long prec);

void fmprb_poly_sin_cos_series_tangent(fmprb_poly_t s, fmprb_poly_t c,
                                    const fmprb_poly_t h, long n, long prec);

void _fmprb_poly_sin_cos_series(fmprb_ptr s, fmprb_ptr c,
                        const fmprb_srcptr h, long hlen, long len, long prec);

void fmprb_poly_sin_cos_series(fmprb_poly_t s, fmprb_poly_t c,
                                    const fmprb_poly_t h, long n, long prec);

void _fmprb_poly_sin_series(fmprb_ptr g, fmprb_srcptr h, long hlen, long n, long prec);

void fmprb_poly_sin_series(fmprb_poly_t g, const fmprb_poly_t h, long n, long prec);

void _fmprb_poly_cos_series(fmprb_ptr g, fmprb_srcptr h, long hlen, long n, long prec);

void fmprb_poly_cos_series(fmprb_poly_t g, const fmprb_poly_t h, long n, long prec);

void _fmprb_poly_tan_series(fmprb_ptr g, fmprb_srcptr h, long hlen, long len, long prec);

void fmprb_poly_tan_series(fmprb_poly_t g, const fmprb_poly_t h, long n, long prec);

void _fmprb_poly_gamma_series(fmprb_ptr res, fmprb_srcptr h, long hlen, long len, long prec);

void fmprb_poly_gamma_series(fmprb_poly_t res, const fmprb_poly_t f, long n, long prec);

void _fmprb_poly_rgamma_series(fmprb_ptr res, fmprb_srcptr h, long hlen, long len, long prec);

void fmprb_poly_rgamma_series(fmprb_poly_t res, const fmprb_poly_t f, long n, long prec);

void _fmprb_poly_lgamma_series(fmprb_ptr res, fmprb_srcptr h, long hlen, long len, long prec);

void fmprb_poly_lgamma_series(fmprb_poly_t res, const fmprb_poly_t f, long n, long prec);

void _fmprb_poly_rising_ui_series(fmprb_ptr res, fmprb_srcptr f, long flen, ulong r, long trunc, long prec);

void fmprb_poly_rising_ui_series(fmprb_poly_t res, const fmprb_poly_t f, ulong r, long trunc, long prec);

void _fmprb_poly_zeta_series(fmprb_ptr res, fmprb_srcptr h, long hlen, const fmprb_t a, int deflate, long len, long prec);

void fmprb_poly_zeta_series(fmprb_poly_t res, const fmprb_poly_t f, const fmprb_t a, int deflate, long n, long prec);

/* Root-finding */

void _fmprb_poly_newton_convergence_factor(fmpr_t convergence_factor,
    fmprb_srcptr poly, long len,
    const fmprb_t convergence_interval, long prec);

int _fmprb_poly_newton_step(fmprb_t xnew, fmprb_srcptr poly, long len,
    const fmprb_t x,
    const fmprb_t convergence_interval,
    const fmpr_t convergence_factor, long prec);

void _fmprb_poly_newton_refine_root(fmprb_t r, fmprb_srcptr poly,
    long len, const fmprb_t start,
    const fmprb_t convergence_interval,
    const fmpr_t convergence_factor,
    long eval_extra_prec,
    long prec);

/* Macros */

/* counts zero bits in the binary representation of e */
static __inline__ int
n_zerobits(mp_limb_t e)
{
    int zeros = 0;

    while (e > 1)
    {
        zeros += !(e & 1);
        e >>= 1;
    }

    return zeros;
}

static __inline__ long
poly_pow_length(long poly_len, ulong exp, long trunc)
{
    mp_limb_t hi, lo;
    umul_ppmm(hi, lo, poly_len - 1, exp);
    add_ssaaaa(hi, lo, hi, lo, 0, 1);
    if (hi != 0 || lo > (mp_limb_t) LONG_MAX)
        return trunc;
    return FLINT_MIN(lo, trunc);
}

#ifndef NEWTON_INIT

#define NEWTON_INIT(from, to) \
    { \
        long __steps[FLINT_BITS], __i, __from, __to; \
        __steps[__i = 0] = __to = (to); \
        __from = (from); \
        while (__to > __from) \
            __steps[++__i] = (__to = (__to + 1) / 2); \

#define NEWTON_BASECASE(bc_to) { long bc_to = __to;

#define NEWTON_END_BASECASE }

#define NEWTON_LOOP(step_from, step_to) \
        { \
            for (__i--; __i >= 0; __i--) \
            { \
                long step_from = __steps[__i+1]; \
                long step_to = __steps[__i]; \

#define NEWTON_END_LOOP }}

#define NEWTON_END }

#endif

#endif

