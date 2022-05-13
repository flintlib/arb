/*
    Copyright (C) 2022 Albin Ahlbäck
    Copyright (C) 2022 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "flint/flint.h"

#ifndef ARB_DEFS_H
#define ARB_DEFS_H

#ifdef __cplusplus
extern "C" {
#endif

#define __ARB_VERSION 2
#define __ARB_VERSION_MINOR 22
#define __ARB_VERSION_PATCHLEVEL 1

#define ARB_VERSION "2.22.1"

#define __ARB_RELEASE (__ARB_VERSION            * 10000 +   \
                       __ARB_VERSION_MINOR      * 100   +   \
                       __ARB_VERSION_PATCHLEVEL)


#define TLS_PREFIX FLINT_TLS_PREFIX

#if defined(_MSC_VER) && defined(ARB_BUILD_DLL)
#define ARB_DLL __declspec(dllexport)
#else
#define ARB_DLL FLINT_DLL
#endif


#ifndef flint_abort
#if __FLINT_RELEASE <= 20502
#define flint_abort abort
#endif
#endif

#if __FLINT_RELEASE < 20600
#define flint_bitcnt_t ulong
#endif


#define LIMB_ONE ((mp_limb_t) 1)
#define LIMB_ONES (-(mp_limb_t) 1)
#define LIMB_TOP (((mp_limb_t) 1) << (FLINT_BITS - 1))
#define MASK_LIMB(n, c) ((n) & (LIMB_ONES << (c)))


/* FLINT_ABS is unsafe for x = WORD_MIN */
#define UI_ABS_SI(x) (((slong)(x) < 0) ? (-(ulong)(x)) : ((ulong)(x)))


#define nn_mul_2x1(r2, r1, r0, a1, a0, b0)                  \
    do {                                                    \
        mp_limb_t t1;                                       \
        umul_ppmm(r1, r0, a0, b0);                          \
        umul_ppmm(r2, t1, a1, b0);                          \
        add_ssaaaa(r2, r1, r2, r1, 0, t1);                  \
    } while (0)

#define nn_mul_2x2(r3, r2, r1, r0, a1, a0, b1, b0)          \
    do {                                                    \
        mp_limb_t t1, t2, t3;                               \
        umul_ppmm(r1, r0, a0, b0);                          \
        umul_ppmm(r2, t1, a1, b0);                          \
        add_ssaaaa(r2, r1, r2, r1, 0, t1);                  \
        umul_ppmm(t1, t2, a0, b1);                          \
        umul_ppmm(r3, t3, a1, b1);                          \
        add_ssaaaa(r3, t1, r3, t1, 0, t3);                  \
        add_ssaaaa(r2, r1, r2, r1, t1, t2);                 \
        r3 += r2 < t1;                                      \
    } while (0)


#ifndef NEWTON_INIT

#define NEWTON_INIT(from, to)                           \
    {                                                   \
        slong __steps[FLINT_BITS], __i, __from, __to;   \
        __steps[__i = 0] = __to = (to);                 \
        __from = (from);                                \
        while (__to > __from)                           \
            __steps[++__i] = (__to = (__to + 1) / 2);   \

#define NEWTON_BASECASE(bc_to) { slong bc_to = __to;

#define NEWTON_END_BASECASE }

#define NEWTON_LOOP(step_from, step_to)                 \
        {                                               \
            for (__i--; __i >= 0; __i--)                \
            {                                           \
                slong step_from = __steps[__i+1];       \
                slong step_to = __steps[__i];           \

#define NEWTON_END_LOOP }}

#define NEWTON_END }

#endif


ARB_DLL extern const char * arb_version;

double arb_test_multiplier(void);

/* counts zero bits in the binary representation of e */
static __inline__
int
n_zerobits(mp_limb_t e)
{
#ifdef __GMP_SHORT_LIMB
    return FLINT_BITS - __builtin_popcount(e) - __builtin_clz(e);
#else
# ifdef _LONG_LONG_LIMB
    return FLINT_BITS - __builtin_popcountll(e) - __builtin_clzll(e);
# else
    return FLINT_BITS - __builtin_popcountl(e) - __builtin_clzl(e);
# endif
#endif
}


#ifdef __cplusplus
}
#endif

#endif
