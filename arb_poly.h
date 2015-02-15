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

#ifndef ARB_POLY_H
#define ARB_POLY_H

#ifdef ARB_POLY_INLINES_C
#define ARB_POLY_INLINE
#else
#define ARB_POLY_INLINE static __inline__
#endif

#include "arb.h"
#include "acb.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    arb_ptr coeffs;
    long length;
    long alloc;
}
arb_poly_struct;

typedef arb_poly_struct arb_poly_t[1];


/* Memory management */

void arb_poly_init(arb_poly_t poly);

void arb_poly_init2(arb_poly_t poly, long len);

void arb_poly_clear(arb_poly_t poly);

void arb_poly_fit_length(arb_poly_t poly, long len);

void _arb_poly_set_length(arb_poly_t poly, long len);

void _arb_poly_normalise(arb_poly_t poly);

ARB_POLY_INLINE void
arb_poly_swap(arb_poly_t poly1, arb_poly_t poly2)
{
    arb_poly_struct t = *poly1;
    *poly1 = *poly2;
    *poly2 = t;
}

void arb_poly_set(arb_poly_t poly, const arb_poly_t src);

void arb_poly_set_round(arb_poly_t poly, const arb_poly_t src, long prec);

/* Basic manipulation */

ARB_POLY_INLINE long arb_poly_length(const arb_poly_t poly)
{
    return poly->length;
}

ARB_POLY_INLINE long arb_poly_degree(const arb_poly_t poly)
{
    return poly->length - 1;
}

ARB_POLY_INLINE void arb_poly_zero(arb_poly_t poly)
{
    poly->length = 0;
}

ARB_POLY_INLINE void
arb_poly_one(arb_poly_t poly)
{
    arb_poly_fit_length(poly, 1);
    arb_one(poly->coeffs);
    _arb_poly_set_length(poly, 1);
}

void arb_poly_set_coeff_si(arb_poly_t poly, long n, long x);

void arb_poly_set_coeff_arb(arb_poly_t poly, long n, const arb_t x);

void arb_poly_get_coeff_arb(arb_t x, const arb_poly_t poly, long n);

#define arb_poly_get_coeff_ptr(poly, n) \
    ((n) < (poly)->length ? (poly)->coeffs + (n) : NULL)

void _arb_poly_reverse(arb_ptr res, arb_srcptr poly, long len, long n);

void _arb_poly_shift_right(arb_ptr res, arb_srcptr poly, long len, long n);

void arb_poly_shift_right(arb_poly_t res, const arb_poly_t poly, long n);

void _arb_poly_shift_left(arb_ptr res, arb_srcptr poly, long len, long n);

void arb_poly_shift_left(arb_poly_t res, const arb_poly_t poly, long n);

ARB_POLY_INLINE void
arb_poly_truncate(arb_poly_t poly, long newlen)
{
    if (poly->length > newlen)
    {
        long i;
        for (i = newlen; i < poly->length; i++)
            arb_zero(poly->coeffs + i);
        poly->length = newlen;
        _arb_poly_normalise(poly);
    }
}

/* Conversions */

void arb_poly_set_fmpz_poly(arb_poly_t poly, const fmpz_poly_t src, long prec);

void arb_poly_set_fmpq_poly(arb_poly_t poly, const fmpq_poly_t src, long prec);

ARB_POLY_INLINE void
arb_poly_set_arb(arb_poly_t poly, const arb_t c)
{
    arb_poly_fit_length(poly, 1);
    arb_set(poly->coeffs, c);
    _arb_poly_set_length(poly, !arb_is_zero(poly->coeffs));
}

void arb_poly_set_si(arb_poly_t poly, long c);

int arb_poly_get_unique_fmpz_poly(fmpz_poly_t res, const arb_poly_t src);

/* Comparisons */

int arb_poly_contains(const arb_poly_t poly1, const arb_poly_t poly2);

int arb_poly_contains_fmpz_poly(const arb_poly_t poly1, const fmpz_poly_t poly2);

int arb_poly_contains_fmpq_poly(const arb_poly_t poly1, const fmpq_poly_t poly2);

int arb_poly_equal(const arb_poly_t A, const arb_poly_t B);

int _arb_poly_overlaps(arb_srcptr poly1, long len1, arb_srcptr poly2, long len2);

int arb_poly_overlaps(const arb_poly_t poly1, const arb_poly_t poly2);

/* Bounds */

void _arb_poly_majorant(arb_ptr res, arb_srcptr vec, long len, long prec);

void arb_poly_majorant(arb_poly_t res, const arb_poly_t poly, long prec);

/* IO */

void arb_poly_printd(const arb_poly_t poly, long digits);

/* Random generation */

void arb_poly_randtest(arb_poly_t poly, flint_rand_t state, long len, long prec, long mag_bits);

/* Arithmetic */

void
_arb_poly_add(arb_ptr res, arb_srcptr poly1, long len1,
    arb_srcptr poly2, long len2, long prec);

void arb_poly_add(arb_poly_t res, const arb_poly_t poly1,
              const arb_poly_t poly2, long prec);

void arb_poly_add_si(arb_poly_t res, const arb_poly_t poly, long c, long prec);

void _arb_poly_sub(arb_ptr res, arb_srcptr poly1, long len1,
    arb_srcptr poly2, long len2, long prec);

void arb_poly_sub(arb_poly_t res, const arb_poly_t poly1,
              const arb_poly_t poly2, long prec);

ARB_POLY_INLINE void
arb_poly_neg(arb_poly_t res, const arb_poly_t poly)
{
    arb_poly_fit_length(res, poly->length);
    _arb_vec_neg(res->coeffs, poly->coeffs, poly->length);
    _arb_poly_set_length(res, poly->length);
}

ARB_POLY_INLINE void
arb_poly_scalar_mul_2exp_si(arb_poly_t res, const arb_poly_t poly, long c)
{
    arb_poly_fit_length(res, poly->length);
    _arb_vec_scalar_mul_2exp_si(res->coeffs, poly->coeffs, poly->length, c);
    _arb_poly_set_length(res, poly->length);
}

void _arb_poly_mullow_ztrunc(arb_ptr res,
    arb_srcptr poly1, long len1,
    arb_srcptr poly2, long len2, long n, long prec);

void arb_poly_mullow_ztrunc(arb_poly_t res, const arb_poly_t poly1,
                                            const arb_poly_t poly2,
                                                long n, long prec);

void _arb_poly_mullow_classical(arb_ptr res,
    arb_srcptr poly1, long len1,
    arb_srcptr poly2, long len2, long n, long prec);

void arb_poly_mullow_classical(arb_poly_t res, const arb_poly_t poly1,
                                            const arb_poly_t poly2,
                                                long n, long prec);

void _arb_poly_mullow_block(arb_ptr C,
    arb_srcptr A, long lenA,
    arb_srcptr B, long lenB, long n, long prec);

void arb_poly_mullow_block(arb_poly_t res, const arb_poly_t poly1,
              const arb_poly_t poly2, long len, long prec);

void _arb_poly_mullow(arb_ptr C,
    arb_srcptr A, long lenA,
    arb_srcptr B, long lenB, long n, long prec);

void arb_poly_mullow(arb_poly_t res, const arb_poly_t poly1,
              const arb_poly_t poly2, long len, long prec);

void _arb_poly_mul(arb_ptr C,
    arb_srcptr A, long lenA,
    arb_srcptr B, long lenB, long prec);

void arb_poly_mul(arb_poly_t res, const arb_poly_t poly1,
              const arb_poly_t poly2, long prec);

ARB_POLY_INLINE void
_arb_poly_mul_monic(arb_ptr res, arb_srcptr poly1, long len1,
    arb_srcptr poly2, long len2, long prec)
{
    if (len1 + len2 - 2 > 0)
        _arb_poly_mullow(res, poly1, len1, poly2, len2, len1 + len2 - 2, prec);
    arb_one(res + len1 + len2 - 2);
}

void _arb_poly_inv_series(arb_ptr Qinv,
    arb_srcptr Q, long Qlen, long len, long prec);

void arb_poly_inv_series(arb_poly_t Qinv, const arb_poly_t Q, long n, long prec);

void  _arb_poly_div_series(arb_ptr Q, arb_srcptr A, long Alen,
    arb_srcptr B, long Blen, long n, long prec);

void arb_poly_div_series(arb_poly_t Q, const arb_poly_t A, const arb_poly_t B, long n, long prec);

void
_arb_poly_div(arb_ptr Q,
    arb_srcptr A, long lenA,
    arb_srcptr B, long lenB, long prec);

void _arb_poly_divrem(arb_ptr Q, arb_ptr R,
    arb_srcptr A, long lenA,
    arb_srcptr B, long lenB, long prec);

void _arb_poly_rem(arb_ptr R,
    arb_srcptr A, long lenA,
    arb_srcptr B, long lenB, long prec);

int arb_poly_divrem(arb_poly_t Q, arb_poly_t R,
                             const arb_poly_t A, const arb_poly_t B, long prec);

void _arb_poly_div_root(arb_ptr Q, arb_t R, arb_srcptr A,
    long len, const arb_t c, long prec);

/* Product trees */

void _arb_poly_product_roots(arb_ptr poly, arb_srcptr xs, long n, long prec);

void arb_poly_product_roots(arb_poly_t poly, arb_srcptr xs, long n, long prec);

arb_ptr * _arb_poly_tree_alloc(long len);

void _arb_poly_tree_free(arb_ptr * tree, long len);

void _arb_poly_tree_build(arb_ptr * tree, arb_srcptr roots, long len, long prec);

/* Composition */

void _arb_poly_compose(arb_ptr res,
    arb_srcptr poly1, long len1,
    arb_srcptr poly2, long len2, long prec);

void arb_poly_compose(arb_poly_t res,
              const arb_poly_t poly1, const arb_poly_t poly2, long prec);

void _arb_poly_compose_horner(arb_ptr res,
    arb_srcptr poly1, long len1,
    arb_srcptr poly2, long len2, long prec);

void arb_poly_compose_horner(arb_poly_t res,
              const arb_poly_t poly1, const arb_poly_t poly2, long prec);

void _arb_poly_compose_divconquer(arb_ptr res,
    arb_srcptr poly1, long len1,
    arb_srcptr poly2, long len2, long prec);

void arb_poly_compose_divconquer(arb_poly_t res,
              const arb_poly_t poly1, const arb_poly_t poly2, long prec);

void _arb_poly_compose_series_horner(arb_ptr res, arb_srcptr poly1, long len1,
                            arb_srcptr poly2, long len2, long n, long prec);

void arb_poly_compose_series_horner(arb_poly_t res,
                    const arb_poly_t poly1,
                    const arb_poly_t poly2, long n, long prec);

void _arb_poly_compose_series(arb_ptr res, arb_srcptr poly1, long len1,
                            arb_srcptr poly2, long len2, long n, long prec);

void arb_poly_compose_series(arb_poly_t res,
                    const arb_poly_t poly1,
                    const arb_poly_t poly2, long n, long prec);

/* Reversion */

void _arb_poly_revert_series_lagrange(arb_ptr Qinv, arb_srcptr Q, long Qlen, long n, long prec);
void arb_poly_revert_series_lagrange(arb_poly_t Qinv, const arb_poly_t Q, long n, long prec);

void _arb_poly_revert_series_newton(arb_ptr Qinv, arb_srcptr Q, long Qlen, long n, long prec);
void arb_poly_revert_series_newton(arb_poly_t Qinv, const arb_poly_t Q, long n, long prec);

void _arb_poly_revert_series_lagrange_fast(arb_ptr Qinv, arb_srcptr Q, long Qlen, long n, long prec);
void arb_poly_revert_series_lagrange_fast(arb_poly_t Qinv, const arb_poly_t Q, long n, long prec);

void _arb_poly_revert_series(arb_ptr Qinv, arb_srcptr Q, long Qlen, long n, long prec);
void arb_poly_revert_series(arb_poly_t Qinv, const arb_poly_t Q, long n, long prec);

/* Evaluation and interpolation */

void _arb_poly_evaluate_horner(arb_t res, arb_srcptr f, long len, const arb_t a, long prec);
void arb_poly_evaluate_horner(arb_t res, const arb_poly_t f, const arb_t a, long prec);

void _arb_poly_evaluate_rectangular(arb_t y, arb_srcptr poly, long len, const arb_t x, long prec);
void arb_poly_evaluate_rectangular(arb_t res, const arb_poly_t f, const arb_t a, long prec);

void _arb_poly_evaluate(arb_t res, arb_srcptr f, long len, const arb_t a, long prec);
void arb_poly_evaluate(arb_t res, const arb_poly_t f, const arb_t a, long prec);

void _arb_poly_evaluate2_horner(arb_t y, arb_t z, arb_srcptr f, long len, const arb_t x, long prec);
void arb_poly_evaluate2_horner(arb_t y, arb_t z, const arb_poly_t f, const arb_t x, long prec);

void _arb_poly_evaluate2_rectangular(arb_t y, arb_t z, arb_srcptr f, long len, const arb_t x, long prec);
void arb_poly_evaluate2_rectangular(arb_t y, arb_t z, const arb_poly_t f, const arb_t x, long prec);

void _arb_poly_evaluate2(arb_t y, arb_t z, arb_srcptr f, long len, const arb_t x, long prec);
void arb_poly_evaluate2(arb_t y, arb_t z, const arb_poly_t f, const arb_t x, long prec);




void _arb_poly_evaluate_vec_iter(arb_ptr ys, arb_srcptr poly, long plen,
    arb_srcptr xs, long n, long prec);

void arb_poly_evaluate_vec_iter(arb_ptr ys,
        const arb_poly_t poly, arb_srcptr xs, long n, long prec);

void _arb_poly_evaluate_vec_fast_precomp(arb_ptr vs, arb_srcptr poly,
    long plen, arb_ptr * tree, long len, long prec);

void _arb_poly_evaluate_vec_fast(arb_ptr ys, arb_srcptr poly, long plen,
    arb_srcptr xs, long n, long prec);

void arb_poly_evaluate_vec_fast(arb_ptr ys,
        const arb_poly_t poly, arb_srcptr xs, long n, long prec);

void _arb_poly_interpolate_newton(arb_ptr poly, arb_srcptr xs,
    arb_srcptr ys, long n, long prec);

void arb_poly_interpolate_newton(arb_poly_t poly,
    arb_srcptr xs, arb_srcptr ys, long n, long prec);

void
_arb_poly_interpolate_barycentric(arb_ptr poly,
    arb_srcptr xs, arb_srcptr ys, long n, long prec);

void arb_poly_interpolate_barycentric(arb_poly_t poly,
    arb_srcptr xs, arb_srcptr ys, long n, long prec);

void _arb_poly_interpolation_weights(arb_ptr w,
    arb_ptr * tree, long len, long prec);

void _arb_poly_interpolate_fast_precomp(arb_ptr poly,
    arb_srcptr ys, arb_ptr * tree, arb_srcptr weights,
    long len, long prec);

void _arb_poly_interpolate_fast(arb_ptr poly,
    arb_srcptr xs, arb_srcptr ys, long len, long prec);

void arb_poly_interpolate_fast(arb_poly_t poly,
        arb_srcptr xs, arb_srcptr ys, long n, long prec);

/* Derivative and integral */

void _arb_poly_derivative(arb_ptr res, arb_srcptr poly, long len, long prec);

void arb_poly_derivative(arb_poly_t res, const arb_poly_t poly, long prec);

void _arb_poly_integral(arb_ptr res, arb_srcptr poly, long len, long prec);

void arb_poly_integral(arb_poly_t res, const arb_poly_t poly, long prec);

/* Transforms */

void arb_poly_borel_transform(arb_poly_t res, const arb_poly_t poly, long prec);

void _arb_poly_borel_transform(arb_ptr res, arb_srcptr poly, long len, long prec);

void arb_poly_inv_borel_transform(arb_poly_t res, const arb_poly_t poly, long prec);

void _arb_poly_inv_borel_transform(arb_ptr res, arb_srcptr poly, long len, long prec);

void _arb_poly_binomial_transform_basecase(arb_ptr b, arb_srcptr a, long alen, long len, long prec);

void arb_poly_binomial_transform_basecase(arb_poly_t b, const arb_poly_t a, long len, long prec);

void _arb_poly_binomial_transform_convolution(arb_ptr b, arb_srcptr a, long alen, long len, long prec);

void arb_poly_binomial_transform_convolution(arb_poly_t b, const arb_poly_t a, long len, long prec);

void _arb_poly_binomial_transform(arb_ptr b, arb_srcptr a, long alen, long len, long prec);

void arb_poly_binomial_transform(arb_poly_t b, const arb_poly_t a, long len, long prec);

/* Special functions */

void _arb_poly_pow_ui_trunc_binexp(arb_ptr res,
    arb_srcptr f, long flen, ulong exp, long len, long prec);

void arb_poly_pow_ui_trunc_binexp(arb_poly_t res,
    const arb_poly_t poly, ulong exp, long len, long prec);

void _arb_poly_pow_ui(arb_ptr res, arb_srcptr f, long flen, ulong exp, long prec);

void arb_poly_pow_ui(arb_poly_t res, const arb_poly_t poly, ulong exp, long prec);

void _arb_poly_pow_series(arb_ptr h,
    arb_srcptr f, long flen,
    arb_srcptr g, long glen, long len, long prec);

void arb_poly_pow_series(arb_poly_t h,
    const arb_poly_t f, const arb_poly_t g, long len, long prec);

void _arb_poly_pow_arb_series(arb_ptr h,
    arb_srcptr f, long flen, const arb_t g, long len, long prec);

void arb_poly_pow_arb_series(arb_poly_t h,
    const arb_poly_t f, const arb_t g, long len, long prec);

void _arb_poly_binomial_pow_arb_series(arb_ptr h, arb_srcptr f, long flen, const arb_t g, long len, long prec);

void _arb_poly_rsqrt_series(arb_ptr g,
    arb_srcptr h, long hlen, long len, long prec);

void arb_poly_rsqrt_series(arb_poly_t g, const arb_poly_t h, long n, long prec);

void _arb_poly_sqrt_series(arb_ptr g,
    arb_srcptr h, long hlen, long len, long prec);

void arb_poly_sqrt_series(arb_poly_t g, const arb_poly_t h, long n, long prec);

void _arb_poly_log_series(arb_ptr res, arb_srcptr f, long flen, long n, long prec);

void arb_poly_log_series(arb_poly_t res, const arb_poly_t f, long n, long prec);

void _arb_poly_atan_series(arb_ptr res, arb_srcptr f, long flen, long n, long prec);

void arb_poly_atan_series(arb_poly_t res, const arb_poly_t f, long n, long prec);

void _arb_poly_asin_series(arb_ptr res, arb_srcptr f, long flen, long n, long prec);

void arb_poly_asin_series(arb_poly_t res, const arb_poly_t f, long n, long prec);

void _arb_poly_acos_series(arb_ptr res, arb_srcptr f, long flen, long n, long prec);

void arb_poly_acos_series(arb_poly_t res, const arb_poly_t f, long n, long prec);

void _arb_poly_exp_series_basecase(arb_ptr f,
        arb_srcptr h, long hlen, long n, long prec);

void arb_poly_exp_series_basecase(arb_poly_t f, const arb_poly_t h, long n, long prec);

void _arb_poly_exp_series(arb_ptr f, arb_srcptr h, long hlen, long n, long prec);

void arb_poly_exp_series(arb_poly_t f, const arb_poly_t h, long n, long prec);

void _arb_poly_sin_cos_series_basecase(arb_ptr s,
                                    arb_ptr c, arb_srcptr h, long hlen, long n, long prec, int times_pi);

void arb_poly_sin_cos_series_basecase(arb_poly_t s, arb_poly_t c,
        const arb_poly_t h, long n, long prec, int times_pi);

void _arb_poly_sin_cos_series_tangent(arb_ptr s, arb_ptr c,
                        const arb_srcptr h, long hlen, long len, long prec, int times_pi);

void arb_poly_sin_cos_series_tangent(arb_poly_t s, arb_poly_t c,
                                    const arb_poly_t h, long n, long prec, int times_pi);

void _arb_poly_sin_cos_series(arb_ptr s, arb_ptr c,
                        const arb_srcptr h, long hlen, long len, long prec);

void arb_poly_sin_cos_series(arb_poly_t s, arb_poly_t c,
                                    const arb_poly_t h, long n, long prec);

void _arb_poly_sin_cos_pi_series(arb_ptr s, arb_ptr c,
                        const arb_srcptr h, long hlen, long len, long prec);

void arb_poly_sin_cos_pi_series(arb_poly_t s, arb_poly_t c,
                                    const arb_poly_t h, long n, long prec);

void _arb_poly_sin_series(arb_ptr g, arb_srcptr h, long hlen, long n, long prec);

void arb_poly_sin_series(arb_poly_t g, const arb_poly_t h, long n, long prec);

void _arb_poly_cos_series(arb_ptr g, arb_srcptr h, long hlen, long n, long prec);

void arb_poly_cos_series(arb_poly_t g, const arb_poly_t h, long n, long prec);

void _arb_poly_sin_pi_series(arb_ptr g, arb_srcptr h, long hlen, long n, long prec);

void arb_poly_sin_pi_series(arb_poly_t g, const arb_poly_t h, long n, long prec);

void _arb_poly_cos_pi_series(arb_ptr g, arb_srcptr h, long hlen, long n, long prec);

void arb_poly_cos_pi_series(arb_poly_t g, const arb_poly_t h, long n, long prec);

void _arb_poly_tan_series(arb_ptr g, arb_srcptr h, long hlen, long len, long prec);

void arb_poly_tan_series(arb_poly_t g, const arb_poly_t h, long n, long prec);

void _arb_poly_compose_series_brent_kung(arb_ptr res, arb_srcptr poly1, long len1,
                            arb_srcptr poly2, long len2, long n, long prec);

void arb_poly_compose_series_brent_kung(arb_poly_t res,
                    const arb_poly_t poly1,
                    const arb_poly_t poly2, long n, long prec);


void _arb_poly_evaluate_acb_horner(acb_t res, arb_srcptr f, long len, const acb_t x, long prec);
void arb_poly_evaluate_acb_horner(acb_t res, const arb_poly_t f, const acb_t a, long prec);

void _arb_poly_evaluate_acb_rectangular(acb_t y, arb_srcptr poly, long len, const acb_t x, long prec);
void arb_poly_evaluate_acb_rectangular(acb_t res, const arb_poly_t f, const acb_t a, long prec);

void _arb_poly_evaluate_acb(acb_t res, arb_srcptr f, long len, const acb_t x, long prec);
void arb_poly_evaluate_acb(acb_t res, const arb_poly_t f, const acb_t a, long prec);

void _arb_poly_evaluate2_acb_horner(acb_t y, acb_t z, arb_srcptr f, long len, const acb_t x, long prec);
void arb_poly_evaluate2_acb_horner(acb_t y, acb_t z, const arb_poly_t f, const acb_t x, long prec);

void _arb_poly_evaluate2_acb_rectangular(acb_t y, acb_t z, arb_srcptr f, long len, const acb_t x, long prec);
void arb_poly_evaluate2_acb_rectangular(acb_t y, acb_t z, const arb_poly_t f, const acb_t x, long prec);

void _arb_poly_evaluate2_acb(acb_t y, acb_t z, arb_srcptr f, long len, const acb_t x, long prec);
void arb_poly_evaluate2_acb(acb_t y, acb_t z, const arb_poly_t f, const acb_t x, long prec);

void _arb_poly_gamma_series(arb_ptr res, arb_srcptr h, long hlen, long len, long prec);
void arb_poly_gamma_series(arb_poly_t res, const arb_poly_t f, long n, long prec);

void _arb_poly_rgamma_series(arb_ptr res, arb_srcptr h, long hlen, long len, long prec);
void arb_poly_rgamma_series(arb_poly_t res, const arb_poly_t f, long n, long prec);

void _arb_poly_lgamma_series(arb_ptr res, arb_srcptr h, long hlen, long len, long prec);
void arb_poly_lgamma_series(arb_poly_t res, const arb_poly_t f, long n, long prec);

void _arb_poly_rising_ui_series(arb_ptr res, arb_srcptr f, long flen, ulong r, long trunc, long prec);
void arb_poly_rising_ui_series(arb_poly_t res, const arb_poly_t f, ulong r, long trunc, long prec);

void _arb_poly_zeta_series(arb_ptr res, arb_srcptr h, long hlen, const arb_t a, int deflate, long len, long prec);
void arb_poly_zeta_series(arb_poly_t res, const arb_poly_t f, const arb_t a, int deflate, long n, long prec);

void _arb_poly_riemann_siegel_theta_series(arb_ptr res, arb_srcptr h, long hlen, long len, long prec);
void arb_poly_riemann_siegel_theta_series(arb_poly_t res, const arb_poly_t h, long n, long prec);

void _arb_poly_riemann_siegel_z_series(arb_ptr res, arb_srcptr h, long hlen, long len, long prec);
void arb_poly_riemann_siegel_z_series(arb_poly_t res, const arb_poly_t h, long n, long prec);

long _arb_poly_swinnerton_dyer_ui_prec(ulong n);
void _arb_poly_swinnerton_dyer_ui(arb_ptr T, ulong n, long trunc, long prec);
void arb_poly_swinnerton_dyer_ui(arb_poly_t poly, ulong n, long prec);

/* Root-finding */

void _arb_poly_newton_convergence_factor(arf_t convergence_factor,
    arb_srcptr poly, long len,
    const arb_t convergence_interval, long prec);

int _arb_poly_newton_step(arb_t xnew, arb_srcptr poly, long len,
    const arb_t x,
    const arb_t convergence_interval,
    const arf_t convergence_factor, long prec);

void _arb_poly_newton_refine_root(arb_t r, arb_srcptr poly,
    long len, const arb_t start,
    const arb_t convergence_interval,
    const arf_t convergence_factor,
    long eval_extra_prec,
    long prec);

/* Macros */


/* counts zero bits in the binary representation of e */
ARB_POLY_INLINE int
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

ARB_POLY_INLINE long
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

#ifdef __cplusplus
}
#endif

#endif

