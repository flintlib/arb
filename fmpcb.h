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

#ifndef FMPCB_H
#define FMPCB_H

#include "fmpr.h"
#include "fmprb.h"

typedef struct
{
    fmprb_struct real;
    fmprb_struct imag;
}
fmpcb_struct;

typedef fmpcb_struct fmpcb_t[1];

#define fmpcb_realref(x) (&(x)->real)
#define fmpcb_imagref(x) (&(x)->imag)


static __inline__ void
fmpcb_init(fmpcb_t x)
{
    fmprb_init(fmpcb_realref(x));
    fmprb_init(fmpcb_imagref(x));
}

static __inline__ void
fmpcb_clear(fmpcb_t x)
{
    fmprb_clear(fmpcb_realref(x));
    fmprb_clear(fmpcb_imagref(x));
}

static __inline__ fmpcb_struct *
_fmpcb_vec_init(long n)
{
    long i;
    fmpcb_struct * v = flint_malloc(sizeof(fmpcb_struct) * n);

    for (i = 0; i < n; i++)
        fmpcb_init(v + i);

    return v;
}

static __inline__ void
_fmpcb_vec_clear(fmpcb_struct * v, long n)
{
    long i;
    for (i = 0; i < n; i++)
        fmpcb_clear(v + i);
    flint_free(v);
}

static __inline__ int
fmpcb_is_zero(const fmpcb_t z)
{
    return fmprb_is_zero(fmpcb_realref(z)) && fmprb_is_zero(fmpcb_imagref(z));
}

static __inline__ int
fmpcb_is_exact(const fmpcb_t z)
{
    return fmprb_is_exact(fmpcb_realref(z)) && fmprb_is_exact(fmpcb_imagref(z));
}


static __inline__ void
fmpcb_zero(fmpcb_t z)
{
    fmprb_zero(fmpcb_realref(z));
    fmprb_zero(fmpcb_imagref(z));
}

static __inline__ void
fmpcb_one(fmpcb_t z)
{
    fmprb_one(fmpcb_realref(z));
    fmprb_zero(fmpcb_imagref(z));
}

static __inline__ void
fmpcb_onei(fmpcb_t z)
{
    fmprb_zero(fmpcb_realref(z));
    fmprb_one(fmpcb_imagref(z));
}

static __inline__ void
fmpcb_set(fmpcb_t z, const fmpcb_t x)
{
    fmprb_set(fmpcb_realref(z), fmpcb_realref(x));
    fmprb_set(fmpcb_imagref(z), fmpcb_imagref(x));
}

static __inline__ void
fmpcb_swap(fmpcb_t z, fmpcb_t x)
{
    fmprb_swap(fmpcb_realref(z), fmpcb_realref(x));
    fmprb_swap(fmpcb_imagref(z), fmpcb_imagref(x));
}

static __inline__ int
fmpcb_overlaps(const fmpcb_t x, const fmpcb_t y)
{
    return fmprb_overlaps(fmpcb_realref(x), fmpcb_realref(y)) &&
            fmprb_overlaps(fmpcb_imagref(x), fmpcb_imagref(y));
}

void
fmpcb_get_abs_ubound_fmpr(fmpr_t u, const fmpcb_t z, long prec)
{
    fmpr_t v;
    fmpr_init(v);

    fmprb_get_abs_ubound_fmpr(u, fmpcb_realref(z), prec);
    fmprb_get_abs_ubound_fmpr(v, fmpcb_realref(z), prec);

    fmpr_mul(u, u, u, prec, FMPR_RND_UP);
    fmpr_mul(v, v, v, prec, FMPR_RND_UP);
    fmpr_add(u, u, v, prec, FMPR_RND_UP);
    fmpr_sqrt(u, u, prec, FMPR_RND_UP);

    fmpr_clear(v);
}

static __inline__ void
fmpcb_add(fmpcb_t z, const fmpcb_t x, const fmpcb_t y, long prec)
{
    fmprb_add(fmpcb_realref(z), fmpcb_realref(x), fmpcb_realref(y), prec);
    fmprb_add(fmpcb_imagref(z), fmpcb_imagref(x), fmpcb_imagref(y), prec);
}

static __inline__ void
fmpcb_sub(fmpcb_t z, const fmpcb_t x, const fmpcb_t y, long prec)
{
    fmprb_sub(fmpcb_realref(z), fmpcb_realref(x), fmpcb_realref(y), prec);
    fmprb_sub(fmpcb_imagref(z), fmpcb_imagref(x), fmpcb_imagref(y), prec);
}

static __inline__ void
fmpcb_neg(fmpcb_t z, const fmpcb_t x)
{
    fmprb_neg(fmpcb_realref(z), fmpcb_realref(x));
    fmprb_neg(fmpcb_imagref(z), fmpcb_imagref(x));
}

static __inline__ void
fmpcb_mul_ui(fmpcb_t z, const fmpcb_t x, ulong y, long prec)
{
    fmprb_mul_ui(fmpcb_realref(z), fmpcb_realref(x), y, prec);
    fmprb_mul_ui(fmpcb_imagref(z), fmpcb_imagref(x), y, prec);
}

static __inline__ void
fmpcb_mul_fmprb(fmpcb_t z, const fmpcb_t x, const fmprb_t y, long prec)
{
    fmprb_mul(fmpcb_realref(z), fmpcb_realref(x), y, prec);
    fmprb_mul(fmpcb_imagref(z), fmpcb_imagref(x), y, prec);
}

static __inline__ void
fmpcb_mul_onei(fmpcb_t z, const fmpcb_t x)
{
    if (z == x)
    {
        fmprb_swap(fmpcb_realref(z), fmpcb_imagref(z));
        fmprb_neg(fmpcb_realref(z), fmpcb_realref(z));
    }
    else
    {
        fmprb_neg(fmpcb_realref(z), fmpcb_imagref(x));
        fmprb_set(fmpcb_imagref(z), fmpcb_realref(x));
    }
}


static __inline__ void
fmpcb_mul(fmpcb_t z, const fmpcb_t x, const fmpcb_t y, long prec)
{
#define a fmpcb_realref(x)
#define b fmpcb_imagref(x)
#define c fmpcb_realref(y)
#define d fmpcb_imagref(y)

/*
    TODO: use these optimizations, but fix aliasing issues

    if (fmprb_is_zero(b))
    {
        fmpcb_mul_fmprb(z, y, a, prec);
    }
    else if (fmprb_is_zero(d))
    {
        fmpcb_mul_fmprb(z, x, c, prec);
    }
    else if (fmprb_is_zero(a))
    {
        fmpcb_mul_fmprb(z, y, b, prec);
        fmpcb_mul_i(z, z);
    }
    else if (fmprb_is_zero(c))
    {
        fmpcb_mul_fmprb(z, x, d, prec);
        fmpcb_mul_i(z, z);
    }
    else
*/
    {
        fmprb_t t, u, v;

        fmprb_init(t);
        fmprb_init(u);
        fmprb_init(v);

        fmprb_add(t, a, b, prec);
        fmprb_add(u, c, d, prec);
        fmprb_mul(v, t, u, prec);

        fmprb_mul(t, a, c, prec);
        fmprb_mul(u, b, d, prec);

        fmprb_sub(fmpcb_realref(z), t, u, prec);
        fmprb_sub(fmpcb_imagref(z), v, t, prec);
        fmprb_sub(fmpcb_imagref(z), fmpcb_imagref(z), u, prec);

        fmprb_clear(t);
        fmprb_clear(u);
        fmprb_clear(v);
    }

#undef a
#undef b
#undef c
#undef d
}

static __inline__ void
fmpcb_inv(fmpcb_t z, const fmpcb_t x, long prec)
{
    fmprb_t t;
    fmprb_init(t);

#define a fmpcb_realref(x)
#define b fmpcb_imagref(x)

    fmprb_mul(t, a, a, prec);
    fmprb_addmul(t, b, b, prec);

    fmprb_div(fmpcb_realref(z), a, t, prec);
    fmprb_div(fmpcb_imagref(z), b, t, prec);

    fmprb_neg(fmpcb_imagref(z), fmpcb_imagref(z));

#undef a
#undef b

    fmprb_clear(t);
}


static __inline__ void
_fmpcb_vec_zero(fmpcb_struct * A, long n)
{
    long i;
    for (i = 0; i < n; i++)
        fmpcb_zero(A + i);
}

static __inline__ void
_fmpcb_vec_set(fmpcb_struct * res, const fmpcb_struct * vec, long len)
{
    long i;
    for (i = 0; i < len; i++)
        fmpcb_set(res + i, vec + i);
}

void fmpcb_printd(const fmpcb_t z, long digits);

#endif
