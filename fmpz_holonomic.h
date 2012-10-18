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

#ifndef FMPZ_HOLONOMIC_H
#define FMPZ_HOLONOMIC_H

#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mat.h"
#include "fmprb.h"
#include "fmprb_mat.h"

static __inline__ void
fmpz_poly_set_si2(fmpz_poly_t poly, long c0, long c1)
{
    fmpz_poly_fit_length(poly, 2);
    fmpz_poly_set_si(poly, c0);
    fmpz_poly_set_coeff_si(poly, 1, c1);
}

static __inline__ void
fmpz_poly_set_si3(fmpz_poly_t poly, long c0, long c1, long c2)
{
    fmpz_poly_fit_length(poly, 3);
    fmpz_poly_set_si(poly, c0);
    fmpz_poly_set_coeff_si(poly, 1, c1);
    fmpz_poly_set_coeff_si(poly, 2, c2);
}

/*******************************************************************************

        Memory management

*******************************************************************************/

typedef struct
{
    fmpz_poly_struct * coeffs;
    long length;
    long alloc;
}
fmpz_holonomic_struct;

typedef fmpz_holonomic_struct fmpz_holonomic_t[1];


void fmpz_holonomic_init(fmpz_holonomic_t op);

void fmpz_holonomic_clear(fmpz_holonomic_t op);

void fmpz_holonomic_fit_length(fmpz_holonomic_t op, long len);

void _fmpz_holonomic_set_length(fmpz_holonomic_t op, long len);


/*******************************************************************************

        Basic operations

*******************************************************************************/

static __inline__ void
fmpz_holonomic_set(fmpz_holonomic_t rop, const fmpz_holonomic_t op)
{
    if (rop != op)
    {
        long i;

        fmpz_holonomic_fit_length(rop, op->length);
        for (i = 0; i < op->length; i++)
            fmpz_poly_set(rop->coeffs + i, op->coeffs + i);
        _fmpz_holonomic_set_length(rop, op->length);
    }
}

static __inline__ int
fmpz_holonomic_equal(const fmpz_holonomic_t op1, const fmpz_holonomic_t op2)
{
    long i;

    if (op1->length != op2->length)
        return 0;

    for (i = 0; i < op1->length; i++)
        if (!fmpz_poly_equal(op1->coeffs + i, op2->coeffs + i))
            return 0;

    return 1;
}

static __inline__ void
fmpz_holonomic_one(fmpz_holonomic_t op)
{
    fmpz_holonomic_fit_length(op, 1);
    fmpz_poly_one(op->coeffs + 0);
    _fmpz_holonomic_set_length(op, 1);
}

static __inline__ long
fmpz_holonomic_order(const fmpz_holonomic_t op)
{
    return op->length - 1;
}

static __inline__ long
fmpz_holonomic_degree(const fmpz_holonomic_t op)
{
    long i, d = -1;

    for (i = 0; i < op->length; i++)
        d = FLINT_MAX(d, fmpz_poly_degree(op->coeffs + i));

    return d;
}

static __inline__ int
fmpz_holonomic_seq_is_constant(const fmpz_holonomic_t op)
{
    return fmpz_holonomic_order(op) == 0;
}

static __inline__ int
fmpz_holonomic_seq_is_cfinite(const fmpz_holonomic_t op)
{
    long i;
    for (i = 0; i < op->length; i++)
        if (fmpz_poly_degree(op->coeffs + i) > 0)
            return 0;
    return 1;
}

static __inline__ int
fmpz_holonomic_seq_is_hypgeom(const fmpz_holonomic_t op)
{
    return fmpz_holonomic_order(op) == 1;
}


/*******************************************************************************

        Input and output

*******************************************************************************/

void fmpz_holonomic_print(const fmpz_holonomic_t op, const char * x, const char * d);


/*******************************************************************************

        Normalisation

*******************************************************************************/

void fmpz_holonomic_normalise_leading(fmpz_holonomic_t op);

void fmpz_holonomic_normalise_sign(fmpz_holonomic_t op);

void fmpz_holonomic_normalise_content(fmpz_holonomic_t op);

void fmpz_holonomic_seq_normalise_trailing(fmpz_holonomic_t op);

static __inline__ void
fmpz_holonomic_seq_normalise(fmpz_holonomic_t op)
{
    fmpz_holonomic_normalise_leading(op);
    fmpz_holonomic_seq_normalise_trailing(op);
    fmpz_holonomic_normalise_content(op);
    fmpz_holonomic_normalise_sign(op);
}

/*******************************************************************************

        Random generation

*******************************************************************************/

void fmpz_holonomic_randtest(fmpz_holonomic_t op, flint_rand_t state, long r, long d, long b);

/*******************************************************************************

        Shifting

*******************************************************************************/

void fmpz_holonomic_shift_fmpz(fmpz_holonomic_t res, const fmpz_holonomic_t op, const fmpz_t s);

void fmpz_holonomic_shift_fmpq(fmpz_holonomic_t res, const fmpz_holonomic_t op, const fmpq_t s);

void fmpz_holonomic_shift_si(fmpz_holonomic_t res, const fmpz_holonomic_t op, long s);

/*******************************************************************************

        Special sequences

*******************************************************************************/

void fmpz_holonomic_seq_set_const(fmpz_holonomic_t op);

void fmpz_holonomic_seq_set_fmpz_pow(fmpz_holonomic_t op, const fmpz_t c);

void fmpz_holonomic_seq_set_fmpq_pow(fmpz_holonomic_t op, const fmpq_t c);

void fmpz_holonomic_seq_set_factorial(fmpz_holonomic_t op);

void fmpz_holonomic_seq_set_harmonic(fmpz_holonomic_t op);

void fmpz_holonomic_seq_set_fibonacci(fmpz_holonomic_t op);


/*******************************************************************************

        Special functions

*******************************************************************************/

void fmpz_holonomic_fun_set_pow_fmpz(fmpz_holonomic_t op, const fmpz_t e);

void fmpz_holonomic_fun_set_pow_fmpq(fmpz_holonomic_t op, const fmpq_t e);

void fmpz_holonomic_fun_set_exp(fmpz_holonomic_t op);

void fmpz_holonomic_fun_set_sin_cos(fmpz_holonomic_t op);

void fmpz_holonomic_fun_set_log(fmpz_holonomic_t op);

void fmpz_holonomic_fun_set_atan(fmpz_holonomic_t op);

void fmpz_holonomic_fun_set_asin_acos(fmpz_holonomic_t op);

void fmpz_holonomic_fun_set_erf(fmpz_holonomic_t op);


/*******************************************************************************

        Sequence transformations

*******************************************************************************/

void fmpz_holonomic_seq_pow_si(fmpz_holonomic_t res, const fmpz_holonomic_t op, long e);

void fmpz_holonomic_seq_mul(fmpz_holonomic_t res, const fmpz_holonomic_t op1, const fmpz_holonomic_t op2);

void fmpz_holonomic_seq_reverse(fmpz_holonomic_t res, const fmpz_holonomic_t op);

void fmpz_holonomic_seq_section(fmpz_holonomic_t res, const fmpz_holonomic_t op, long m);


/*******************************************************************************

        Taylor series

*******************************************************************************/

void fmpz_holonomic_get_series(fmpz_holonomic_t re, const fmpz_holonomic_t de);


/*******************************************************************************

        Sequence evaluation

*******************************************************************************/

void fmpz_holonomic_forward_fmpz_mat(fmpz_mat_t M, fmpz_t Q, const fmpz_holonomic_t op, long start, long n);

void fmpz_holonomic_get_nth_fmpz(fmpz_t res, const fmpz_holonomic_t op, const fmpz * initial, long n0, long n);

void fmpz_holonomic_get_nth_fmpq(fmpq_t res, const fmpz_holonomic_t op, const fmpq * initial, long n0, long n);

void fmpz_holonomic_forward_fmprb_mat(fmprb_mat_t M, fmprb_t Q, const fmpz_holonomic_t op, long start, long n, long prec);

#endif

