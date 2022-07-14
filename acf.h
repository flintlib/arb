/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#ifndef ACF_H
#define ACF_H

#ifdef ACF_INLINES_C
#define ACF_INLINE
#else
#define ACF_INLINE static __inline__
#endif

#include "arf.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    arf_struct real;
    arf_struct imag;
}
acf_struct;

typedef acf_struct acf_t[1];
typedef acf_struct * acf_ptr;
typedef const acf_struct * acf_srcptr;

#define acf_realref(x) (&(x)->real)
#define acf_imagref(x) (&(x)->imag)

ACF_INLINE void
acf_init(acf_t x)
{
    arf_init(acf_realref(x));
    arf_init(acf_imagref(x));
}

ACF_INLINE void
acf_clear(acf_t x)
{
    arf_clear(acf_realref(x));
    arf_clear(acf_imagref(x));
}

ACF_INLINE arf_ptr acf_real_ptr(acf_t z) { return acf_realref(z); }
ACF_INLINE arf_ptr acf_imag_ptr(acf_t z) { return acf_imagref(z); }

ACF_INLINE void
acf_set(acf_t z, const acf_t x)
{
    arf_set(acf_realref(z), acf_realref(x));
    arf_set(acf_imagref(z), acf_imagref(x));
}

ACF_INLINE void
acf_swap(acf_t z, acf_t x)
{
    arf_swap(acf_realref(z), acf_realref(x));
    arf_swap(acf_imagref(z), acf_imagref(x));
}

ACF_INLINE int
acf_equal(const acf_t x, const acf_t y)
{
    return arf_equal(acf_realref(x), acf_realref(y)) &&
           arf_equal(acf_imagref(x), acf_imagref(y));
}

ACF_INLINE slong
acf_allocated_bytes(const acf_t x)
{
    return arf_allocated_bytes(acf_realref(x)) + arf_allocated_bytes(acf_imagref(x));
}

#ifdef __cplusplus
}
#endif

#endif
