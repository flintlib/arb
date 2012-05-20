/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include <alloca.h>
#include "mpr.h"

#define EXP_CACHE_PREC_LIMBS 6
#define EXP_CACHE1_BITS 8
#define EXP_CACHE2_BITS 8


mp_limb_t exp_cache_1[1 << EXP_CACHE1_BITS][EXP_CACHE_PREC_LIMBS + 1];
mp_limb_t exp_cache_2[1 << EXP_CACHE1_BITS][EXP_CACHE_PREC_LIMBS + 1];

int exp_cache_initialised = 0;

mp_size_t
fixed_set_mpfr(mp_ptr y, mpfr_srcptr t, mp_size_t prec)
{
    mp_size_t n;
    mpz_t z;
    mpfr_t u;

    mpz_init(z);
    mpfr_init2(u, mpfr_get_prec(t));

    mpfr_mul_2exp(u, t, prec * FLINT_BITS, GMP_RNDD);
    mpfr_get_z(z, u, GMP_RNDD);
    n = z->_mp_size;

    mpn_copyi(y, z->_mp_d, n);
    mpn_zero(y + n, prec - n);

    mpz_clear(z);
    mpfr_clear(u);

    return n;
}

void
exp_cache_init()
{
    int i, prec;
    mpfr_t h, exph;

    prec = EXP_CACHE_PREC_LIMBS * FLINT_BITS + 1;

    mpfr_init2(h, prec + 64);
    mpfr_init2(exph, prec + 64);
    mpfr_set_ui_2exp(h, 1, -EXP_CACHE1_BITS, GMP_RNDN);
    mpfr_exp(exph, h, GMP_RNDN);
    mpfr_set_ui(h, 1, GMP_RNDN);


    mpfr_mul_2exp(h, h, 64 + 16, GMP_RNDN);
    mpfr_div_ui(h, h, 2432902008176640000UL, GMP_RNDN);

    for (i = 0; i < (1 << EXP_CACHE1_BITS); i++)
    {
        fixed_set_mpfr(exp_cache_1[i], h, EXP_CACHE_PREC_LIMBS);
        mpfr_mul(h, h, exph, GMP_RNDN);
    }

    mpfr_set_ui_2exp(h, 1, -EXP_CACHE1_BITS-EXP_CACHE2_BITS, GMP_RNDN);
    mpfr_exp(exph, h, GMP_RNDN);
    mpfr_set_ui(h, 1, GMP_RNDN);



    mpfr_mul_2exp(h, h, 16, GMP_RNDN);

    for (i = 0; i < (1 << EXP_CACHE2_BITS); i++)
    {
        fixed_set_mpfr(exp_cache_2[i], h, EXP_CACHE_PREC_LIMBS);
        mpfr_mul(h, h, exph, GMP_RNDN);
    }

    mpfr_clear(h);
    mpfr_clear(exph);

    exp_cache_initialised = 1;
}



static const mp_limb_t exp_coeffs[] =
{
#if FLINT64
    2432902008176640000UL,
    2432902008176640000UL,
    1216451004088320000UL,
    405483668029440000UL,
    101370917007360000UL,
    20274183401472000UL,
    3379030566912000UL,
    482718652416000UL,
    60339831552000UL,
    6704425728000UL,
    670442572800UL,
    60949324800UL,
    5079110400UL,
    390700800UL,
    27907200UL,
    1860480UL,
    116280UL,
    6840UL,
    380UL,
    20UL,
    1UL,
#else
    0UL,
    479001600UL,
    239500800UL,
    79833600UL,
    19958400UL,
    3991680UL,
    665280UL,
    95040UL,
    11880UL,
    1320UL,
    132UL,
    12UL,
    1UL,
#endif
};

#define EXP_DIVFREE_MAXTERMS (sizeof(exp_coeffs) / sizeof(mp_limb_t))


static const unsigned char fac_bits[] =
{
    1, 1, 2, 3, 5, 7, 10, 13, 16, 19, 22, 26, 29, 33, 37, 41, 45, 49,
    53, 57, 62, 66, 70, 75, 80, 84, 89, 94, 98, 103, 108, 113, 118, 123,
    128, 133, 139, 144, 149, 154, 160, 165, 170, 176, 181, 187, 192, 198,
    203, 209, 215, 220, 226, 232, 238, 243, 249, 255
};


static __inline__ long
exp_needed_terms(long reduced, long tol_bits)
{
    int i;

    i = 2 + (tol_bits + 1) / (1 + reduced);
    while (reduced*i + fac_bits[i] - 1 > tol_bits) i--;

    return i + 1;
}

void
mpr_exp_basecase(mp_ptr y, mp_srcptr x, long limbs, mp_bitcnt_t tol_bits)
{
    mp_limb_t top;
    long terms;
    int i1, i2;

    if (exp_cache_initialised == 0)
        exp_cache_init();

    top = x[limbs - 1];
    i1 = top >> (FLINT_BITS - EXP_CACHE1_BITS);
    i2 = (top >> (FLINT_BITS - (EXP_CACHE1_BITS + EXP_CACHE2_BITS))) & \
        ((1UL<<EXP_CACHE2_BITS) - 1);

    terms = exp_needed_terms(EXP_CACHE1_BITS + EXP_CACHE2_BITS, tol_bits);

    if (terms <= EXP_DIVFREE_MAXTERMS)
    {
        mp_limb_t t[EXP_CACHE_PREC_LIMBS*2 + 2];
        mp_limb_t u[EXP_CACHE_PREC_LIMBS*2 + 2];

        mpn_copyi(t, x, limbs - 1);
        t[limbs - 1] = (top << (EXP_CACHE1_BITS + EXP_CACHE2_BITS)) >> \
                            (EXP_CACHE1_BITS + EXP_CACHE2_BITS);

        mpr_polyval_1(u, t, limbs, exp_coeffs, terms);
        mpn_rshift(u, u, limbs + 1, 48);

        /* exp(x1+x2+t) = exp(x1)*exp(x2)*exp(t) */
        mpn_mul_n(t, u,         exp_cache_1[i1] + (EXP_CACHE_PREC_LIMBS - limbs), limbs + 1);
        mpn_mul_n(u, t + limbs, exp_cache_2[i2] + (EXP_CACHE_PREC_LIMBS - limbs), limbs + 1);

        /* divide out inflation
        mpn_rshift(u + limbs, u + limbs, limbs + 1, 48);
        mpn_copyi(y, u + limbs, limbs + 1);
         */

        mpn_rshift(y, u + limbs, limbs + 1, 48);
        return;
    }

    printf("exp: not implemented: %ld terms\n", terms);
    abort();
}
