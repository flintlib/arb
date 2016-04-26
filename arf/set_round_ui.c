/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arf.h"

/* Top-aligns the nonzero limb v and rounds it to prec bit. */
#define ARF_NORMALISE_ROUND_LIMB(inexact, exp, v,               \
                                    sgnbit, prec, rnd)          \
    do {                                                        \
        count_leading_zeros(exp, v);                            \
        v <<= exp;                                              \
        exp = FLINT_BITS - exp;                                 \
        if (prec >= exp)                                        \
        {                                                       \
            inexact = 0;                                        \
        }                                                       \
        else                                                    \
        {                                                       \
            mp_limb_t __t = v;                                  \
            v = MASK_LIMB(v, FLINT_BITS - prec);                \
            inexact = (__t != v);                               \
            if (inexact && arf_rounds_up(rnd, sgnbit))          \
            {                                                   \
                v += (LIMB_ONE << (FLINT_BITS - prec));         \
                if (v == 0)                                     \
                {                                               \
                    v = LIMB_TOP;                               \
                    exp++;                                      \
                }                                               \
            }                                                   \
        }                                                       \
    }                                                           \
    while (0)

int
_arf_set_round_ui(arf_t x, ulong v, int sgnbit, slong prec, arf_rnd_t rnd)
{
    if (v == 0)
    {
        arf_zero(x);
        return 0;
    }
    else
    {
        int exp, inexact;

        ARF_NORMALISE_ROUND_LIMB(inexact, exp, v, sgnbit, prec, rnd);

        _fmpz_demote(ARF_EXPREF(x));
        ARF_EXP(x) = exp;
        ARF_DEMOTE(x);
        ARF_XSIZE(x) = ARF_MAKE_XSIZE(1, sgnbit);
        ARF_NOPTR_D(x)[0] = v;
        return inexact;
    }
}

