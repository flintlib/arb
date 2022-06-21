/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef FMPZI_H
#define FMPZI_H

#ifdef FMPZI_INLINES_C
#define FMPZI_INLINE
#else
#define FMPZI_INLINE static __inline__
#endif

#include "flint/fmpz.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    fmpz a;
    fmpz b;
}
fmpzi_struct;

typedef fmpzi_struct fmpzi_t[1];

#define fmpzi_realref(x) (&((x)->a))
#define fmpzi_imagref(x) (&((x)->b))

FMPZI_INLINE void
fmpzi_init(fmpzi_t x)
{
    fmpz_init(fmpzi_realref(x));
    fmpz_init(fmpzi_imagref(x));
}

FMPZI_INLINE void
fmpzi_clear(fmpzi_t x)
{
    fmpz_clear(fmpzi_realref(x));
    fmpz_clear(fmpzi_imagref(x));
}

FMPZI_INLINE int
fmpzi_equal(const fmpzi_t x, const fmpzi_t y)
{
    return fmpz_equal(fmpzi_realref(x), fmpzi_realref(y)) &&
           fmpz_equal(fmpzi_imagref(x), fmpzi_imagref(y));
}

FMPZI_INLINE void
fmpzi_zero(fmpzi_t x)
{
    fmpz_zero(fmpzi_realref(x));
    fmpz_zero(fmpzi_imagref(x));
}

FMPZI_INLINE void
fmpzi_one(fmpzi_t x)
{
    fmpz_one(fmpzi_realref(x));
    fmpz_zero(fmpzi_imagref(x));
}

FMPZI_INLINE void
fmpzi_set(fmpzi_t res, const fmpzi_t x)
{
    fmpz_set(fmpzi_realref(res), fmpzi_realref(x));
    fmpz_set(fmpzi_imagref(res), fmpzi_imagref(x));
}

FMPZI_INLINE void
fmpzi_conj(fmpzi_t res, const fmpzi_t x)
{
    fmpz_set(fmpzi_realref(res), fmpzi_realref(x));
    fmpz_neg(fmpzi_imagref(res), fmpzi_imagref(x));
}

FMPZI_INLINE void
fmpzi_swap(fmpzi_t x, fmpzi_t y)
{
    fmpzi_struct t;
    t = *x;
    *x = *y;
    *y = t;
}

FMPZI_INLINE void
fmpzi_print(const fmpzi_t x)
{
    fmpz_print(fmpzi_realref(x));
    if (fmpz_sgn(fmpzi_imagref(x)) >= 0)
        flint_printf("+");
    fmpz_print(fmpzi_imagref(x));
    flint_printf("*I");
}

FMPZI_INLINE void
fmpzi_set_si_si(fmpzi_t res, slong a, slong b)
{
    fmpz_set_si(fmpzi_realref(res), a);
    fmpz_set_si(fmpzi_imagref(res), b);
}

FMPZI_INLINE void
fmpzi_randtest(fmpzi_t res, flint_rand_t state, mp_bitcnt_t bits)
{
    fmpz_randtest(fmpzi_realref(res), state, bits);
    fmpz_randtest(fmpzi_imagref(res), state, bits);
}

FMPZI_INLINE void
fmpzi_neg(fmpzi_t res, const fmpzi_t x)
{
    fmpz_neg(fmpzi_realref(res), fmpzi_realref(x));
    fmpz_neg(fmpzi_imagref(res), fmpzi_imagref(x));
}

FMPZI_INLINE void
fmpzi_add(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)
{
    fmpz_add(fmpzi_realref(res), fmpzi_realref(x), fmpzi_realref(y));
    fmpz_add(fmpzi_imagref(res), fmpzi_imagref(x), fmpzi_imagref(y));
}

FMPZI_INLINE void
fmpzi_sub(fmpzi_t res, const fmpzi_t x, const fmpzi_t y)
{
    fmpz_sub(fmpzi_realref(res), fmpzi_realref(x), fmpzi_realref(y));
    fmpz_sub(fmpzi_imagref(res), fmpzi_imagref(x), fmpzi_imagref(y));
}

void fmpzi_sqr(fmpzi_t res, const fmpzi_t x);
void fmpzi_mul(fmpzi_t res, const fmpzi_t x, const fmpzi_t y);
void fmpzi_pow_ui(fmpzi_t res, const fmpzi_t x, ulong exp);

#ifdef __cplusplus
}
#endif

#endif
