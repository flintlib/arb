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

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include "arf.h"
#include "flint/mpn_extras.h"

#if !defined(__MPIR_VERSION)
#define USE_GMP_DIV_Q 1
void __gmpn_div_q(mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t, mp_ptr);
#else
#define USE_GMP_DIV_Q 0
#endif

void
arf_div_special(arf_t z, const arf_t x, const arf_t y)
{
    if ((arf_is_zero(x) && !arf_is_zero(y) && !arf_is_nan(y)) ||
        (arf_is_inf(y) && !arf_is_special(x)))
    {
        arf_zero(z);
    }
    else if (arf_is_zero(y) || (arf_is_special(x) && arf_is_special(y)) ||
        arf_is_nan(x) || arf_is_nan(y))
    {
        arf_nan(z);
    }
    else if (arf_sgn(x) == arf_sgn(y))
        arf_pos_inf(z);
    else
        arf_neg_inf(z);
}

int
arf_div(arf_ptr z, arf_srcptr x, arf_srcptr y, slong prec, arf_rnd_t rnd)
{
    mp_size_t xn, yn, zn, sn, tn, alloc;
    mp_srcptr xptr, yptr;
    mp_ptr tmp;
    mp_ptr tptr, zptr;
    int inexact;
    slong fix, fix2;
    ARF_MUL_TMP_DECL

    if (arf_is_special(x) || arf_is_special(y))
    {
        arf_div_special(z, x, y);
        return 0;
    }

    ARF_GET_MPN_READONLY(xptr, xn, x);
    ARF_GET_MPN_READONLY(yptr, yn, y);

    /* Division by a power of two */
    if (yn == 1 && yptr[0] == LIMB_TOP)
    {
        fmpz_t t;
        fmpz_init_set(t, ARF_EXPREF(y));

        if (ARF_SGNBIT(y))
            inexact = arf_neg_round(z, x, prec, rnd);
        else
            inexact = arf_set_round(z, x, prec, rnd);

        _fmpz_sub2_fast(ARF_EXPREF(z), ARF_EXPREF(z), t, 1);
        fmpz_clear(t);
        return inexact;
    }

    sn = FLINT_MAX(prec - xn * FLINT_BITS + yn * FLINT_BITS, 0) + 32;
    sn = (sn + FLINT_BITS - 1) / FLINT_BITS;

    /* Need space for numerator (tn), quotient (zn), and one extra limb at the
       top of tptr for the multiplication to check the remainder. */
    tn = xn + sn;
    zn = tn - yn + 1;
    alloc = zn + (tn + 1);
/* need tn + 1 extra temporary limbs, which we store at the end of tptr */
#if USE_GMP_DIV_Q
    alloc += tn + 1;
#endif

    ARF_MUL_TMP_ALLOC(tmp, alloc)

    zptr = tmp;
    tptr = tmp + zn;

    flint_mpn_zero(tptr, sn);
    flint_mpn_copyi(tptr + sn, xptr, xn);
/* uses tn + 1 extra temporary limbs, tn limbs after tptr */
#if USE_GMP_DIV_Q
    __gmpn_div_q(zptr, tptr, tn, yptr, yn, tptr + tn);
#else
    mpn_tdiv_q(zptr, tptr, tn, yptr, yn);
#endif

    if (zptr[zn - 1] == 0)
    {
        zn--;
        fix2 = 0;
    }
    else
    {
        fix2 = FLINT_BITS;
    }

    /* The remainder can only be zero if the last several bits of the
       extended quotient are zero. */
    if ((zptr[0] & ((LIMB_ONE << 24) - 1)) == 0)
    {
        /* The quotient is exact iff multiplying by y gives back the
           input. Note: the multiplication may write sn + xn + 1 limbs, but
           tptr[sn + xn] is guaranteed to be zero in that case since the
           approximate quotient cannot be larger than the true quotient. */
        if (zn >= yn)
            mpn_mul(tptr, zptr, zn, yptr, yn);
        else
            mpn_mul(tptr, yptr, yn, zptr, zn);

        /* The quotient is not exact. Perturbing the approximate quotient
           and rounding gives the correct the result. */
        if (!flint_mpn_zero_p(tptr, sn) ||
            mpn_cmp(tptr + sn, xptr, xn) != 0)
        {
            zptr[0]++;
        }
    }

    inexact = _arf_set_round_mpn(z, &fix, zptr, zn,
            ARF_SGNBIT(x) ^ ARF_SGNBIT(y), prec, rnd);
    _fmpz_sub2_fast(ARF_EXPREF(z), ARF_EXPREF(x), ARF_EXPREF(y), fix + fix2);
    ARF_MUL_TMP_FREE(tmp, alloc)

    return inexact;
}

