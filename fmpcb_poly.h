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

#ifndef FMPCB_POLY_H
#define FMPCB_POLY_H

#include "fmpcb.h"
#include "fmprb_poly.h"

typedef struct
{
    fmpcb_ptr coeffs;
    long length;
    long alloc;
}
fmpcb_poly_struct;

typedef fmpcb_poly_struct fmpcb_poly_t[1];


/* Memory management */

void fmpcb_poly_init(fmpcb_poly_t poly);

void fmpcb_poly_init2(fmpcb_poly_t poly, long len);

void fmpcb_poly_clear(fmpcb_poly_t poly);

void fmpcb_poly_fit_length(fmpcb_poly_t poly, long len);

void _fmpcb_poly_set_length(fmpcb_poly_t poly, long len);

void _fmpcb_poly_normalise(fmpcb_poly_t poly);

static __inline__ void
fmpcb_poly_swap(fmpcb_poly_t poly1, fmpcb_poly_t poly2)
{
    fmpcb_poly_struct t = *poly1;
    *poly1 = *poly2;
    *poly2 = t;
}

static __inline__ long fmpcb_poly_length(const fmpcb_poly_t poly)
{
    return poly->length;
}

static __inline__ void fmpcb_poly_zero(fmpcb_poly_t poly)
{
    poly->length = 0;
}

static __inline__ void
fmpcb_poly_one(fmpcb_poly_t poly)
{
    fmpcb_poly_fit_length(poly, 1);
    fmpcb_one(poly->coeffs);
    _fmpcb_poly_set_length(poly, 1);
}

void fmpcb_poly_printd(const fmpcb_poly_t poly, long digits);

void _fmpcb_poly_evaluate(fmpcb_t res, fmpcb_srcptr f, long len, const fmpcb_t a, long prec);

void fmpcb_poly_evaluate(fmpcb_t res, const fmpcb_poly_t f, const fmpcb_t a, long prec);

void _fmpcb_poly_derivative(fmpcb_ptr res, fmpcb_srcptr poly, long len, long prec);

void fmpcb_poly_derivative(fmpcb_poly_t res, const fmpcb_poly_t poly, long prec);

void fmpcb_poly_set(fmpcb_poly_t dest, const fmpcb_poly_t src);

void fmpcb_poly_set_fmprb_poly(fmpcb_poly_t poly, const fmprb_poly_t re);

void fmpcb_poly_set2_fmprb_poly(fmpcb_poly_t poly, const fmprb_poly_t re, const fmprb_poly_t im);

void fmpcb_poly_set_fmpq_poly(fmpcb_poly_t poly, const fmpq_poly_t re, long prec);

void fmpcb_poly_set2_fmpq_poly(fmpcb_poly_t poly, const fmpq_poly_t re, const fmpq_poly_t im, long prec);

void fmpcb_poly_randtest(fmpcb_poly_t poly, flint_rand_t state, long len, long prec, long mag_bits);

int fmpcb_poly_equal(const fmpcb_poly_t A, const fmpcb_poly_t B);

int fmpcb_poly_contains_fmpq_poly(const fmpcb_poly_t poly1, const fmpq_poly_t poly2);

int _fmpcb_poly_overlaps(fmpcb_srcptr poly1, long len1,
        fmpcb_srcptr poly2, long len2);

int fmpcb_poly_overlaps(const fmpcb_poly_t poly1, const fmpcb_poly_t poly2);

int fmpcb_poly_contains(const fmpcb_poly_t poly1, const fmpcb_poly_t poly2);

void _fmpcb_poly_add(fmpcb_ptr res, fmpcb_srcptr poly1, long len1,
    fmpcb_srcptr poly2, long len2, long prec);

void fmpcb_poly_add(fmpcb_poly_t res, const fmpcb_poly_t poly1,
              const fmpcb_poly_t poly2, long prec);

void _fmpcb_poly_sub(fmpcb_ptr res, fmpcb_srcptr poly1, long len1,
    fmpcb_srcptr poly2, long len2, long prec);

void fmpcb_poly_sub(fmpcb_poly_t res, const fmpcb_poly_t poly1,
              const fmpcb_poly_t poly2, long prec);

void fmpcb_poly_mullow_classical(fmpcb_poly_t res, const fmpcb_poly_t poly1,
                                            const fmpcb_poly_t poly2,
                                                long n, long prec);

void _fmpcb_poly_mullow_classical(fmpcb_ptr res,
    fmpcb_srcptr poly1, long len1,
    fmpcb_srcptr poly2, long len2, long n, long prec);

void _fmpcb_poly_mullow_transpose(fmpcb_ptr res,
    fmpcb_srcptr poly1, long len1,
    fmpcb_srcptr poly2, long len2, long n, long prec);

void fmpcb_poly_mullow_transpose(fmpcb_poly_t res, const fmpcb_poly_t poly1,
                                            const fmpcb_poly_t poly2,
                                                long n, long prec);

void _fmpcb_poly_mullow_transpose_gauss(fmpcb_ptr res,
    fmpcb_srcptr poly1, long len1,
    fmpcb_srcptr poly2, long len2, long n, long prec);

void fmpcb_poly_mullow_transpose_gauss(fmpcb_poly_t res, const fmpcb_poly_t poly1,
                                            const fmpcb_poly_t poly2,
                                                long n, long prec);

void _fmpcb_poly_mullow(fmpcb_ptr res,
    fmpcb_srcptr poly1, long len1,
    fmpcb_srcptr poly2, long len2, long n, long prec);

void fmpcb_poly_mullow(fmpcb_poly_t res, const fmpcb_poly_t poly1,
                                            const fmpcb_poly_t poly2,
                                                long n, long prec);

void _fmpcb_poly_mul(fmpcb_ptr C,
    fmpcb_srcptr A, long lenA,
    fmpcb_srcptr B, long lenB, long prec);

void fmpcb_poly_mul(fmpcb_poly_t res, const fmpcb_poly_t poly1,
              const fmpcb_poly_t poly2, long prec);

static __inline__ void
_fmpcb_poly_mul_monic(fmpcb_ptr res, fmpcb_srcptr poly1, long len1,
    fmpcb_srcptr poly2, long len2, long prec)
{
    if (len1 + len2 - 2 > 0)
        _fmpcb_poly_mullow(res, poly1, len1, poly2, len2, len1 + len2 - 2, prec);
    fmpcb_one(res + len1 + len2 - 2);
}

void _fmpcb_poly_inv_series(fmpcb_ptr Qinv, fmpcb_srcptr Q, long Qlen, long len, long prec);

void fmpcb_poly_inv_series(fmpcb_poly_t Qinv, const fmpcb_poly_t Q, long n, long prec);

void  _fmpcb_poly_div_series(fmpcb_ptr Q, fmpcb_srcptr A, long Alen,
    fmpcb_srcptr B, long Blen, long n, long prec);

void fmpcb_poly_div_series(fmpcb_poly_t Q, const fmpcb_poly_t A, const fmpcb_poly_t B, long n, long prec);

void _fmpcb_poly_reverse(fmpcb_ptr res, fmpcb_srcptr poly, long len, long n);

void _fmpcb_poly_div(fmpcb_ptr Q,
    fmpcb_srcptr A, long lenA,
    fmpcb_srcptr B, long lenB, long prec);

void _fmpcb_poly_divrem(fmpcb_ptr Q, fmpcb_ptr R,
    fmpcb_srcptr A, long lenA,
    fmpcb_srcptr B, long lenB, long prec);

void _fmpcb_poly_rem(fmpcb_ptr R,
    fmpcb_srcptr A, long lenA,
    fmpcb_srcptr B, long lenB, long prec);

void fmpcb_poly_divrem(fmpcb_poly_t Q, fmpcb_poly_t R,
                             const fmpcb_poly_t A, const fmpcb_poly_t B, long prec);

void _fmpcb_poly_div_root(fmpcb_ptr Q, fmpcb_t R, fmpcb_srcptr A,
    long len, const fmpcb_t c, long prec);

void
_fmpcb_poly_evaluate_vec_fast_precomp(fmpcb_ptr vs, fmpcb_srcptr poly,
    long plen, fmpcb_ptr * tree, long len, long prec);

void _fmpcb_poly_evaluate_vec_fast(fmpcb_ptr ys, fmpcb_srcptr poly, long plen,
    fmpcb_srcptr xs, long n, long prec);

void
fmpcb_poly_evaluate_vec_fast(fmpcb_ptr ys,
        const fmpcb_poly_t poly, fmpcb_srcptr xs, long n, long prec);

void
_fmpcb_poly_evaluate_vec_iter(fmpcb_ptr ys, fmpcb_srcptr poly, long plen,
    fmpcb_srcptr xs, long n, long prec);

void
fmpcb_poly_evaluate_vec_iter(fmpcb_ptr ys,
        const fmpcb_poly_t poly, fmpcb_srcptr xs, long n, long prec);

void
_fmpcb_poly_interpolate_barycentric(fmpcb_ptr poly,
    fmpcb_srcptr xs, fmpcb_srcptr ys, long n, long prec);

void
fmpcb_poly_interpolate_barycentric(fmpcb_poly_t poly,
    fmpcb_srcptr xs, fmpcb_srcptr ys, long n, long prec);

void
_fmpcb_poly_interpolation_weights(fmpcb_ptr w,
    fmpcb_ptr * tree, long len, long prec);

void
_fmpcb_poly_interpolate_fast_precomp(fmpcb_ptr poly,
    fmpcb_srcptr ys, fmpcb_ptr * tree, fmpcb_srcptr weights,
    long len, long prec);

void
_fmpcb_poly_interpolate_fast(fmpcb_ptr poly,
    fmpcb_srcptr xs, fmpcb_srcptr ys, long len, long prec);

void
fmpcb_poly_interpolate_fast(fmpcb_poly_t poly,
        fmpcb_srcptr xs, fmpcb_srcptr ys, long n, long prec);

void
_fmpcb_poly_interpolate_newton(fmpcb_ptr poly, fmpcb_srcptr xs,
    fmpcb_srcptr ys, long n, long prec);

void
fmpcb_poly_interpolate_newton(fmpcb_poly_t poly,
    fmpcb_srcptr xs, fmpcb_srcptr ys, long n, long prec);

void
_fmpcb_poly_product_roots(fmpcb_ptr poly, fmpcb_srcptr xs, long n, long prec);

void
fmpcb_poly_product_roots(fmpcb_poly_t poly, fmpcb_ptr xs, long n, long prec);

fmpcb_ptr * _fmpcb_poly_tree_alloc(long len);

void _fmpcb_poly_tree_free(fmpcb_ptr * tree, long len);

void
_fmpcb_poly_tree_build(fmpcb_ptr * tree, fmpcb_srcptr roots, long len, long prec);


void _fmpcb_poly_root_inclusion(fmpcb_t r, const fmpcb_t m,
    fmpcb_srcptr poly,
    fmpcb_srcptr polyder, long len, long prec);

long _fmpcb_poly_validate_roots(fmpcb_ptr roots,
        fmpcb_srcptr poly, long len, long prec);

void _fmpcb_poly_refine_roots_durand_kerner(fmpcb_ptr roots,
        fmpcb_srcptr poly, long len, long prec);

long _fmpcb_poly_find_roots(fmpcb_ptr roots,
    fmpcb_srcptr poly,
    fmpcb_srcptr initial, long len, long maxiter, long prec);

long fmpcb_poly_find_roots(fmpcb_ptr roots,
    const fmpcb_poly_t poly, fmpcb_srcptr initial,
    long maxiter, long prec);

#endif
