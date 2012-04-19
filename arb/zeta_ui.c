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

#include "arb.h"
#include "arith.h"

static void
mpfr_zeta_ui_bernoulli(mpfr_t x, ulong n)
{
    fmpq_t b;
    mpfr_t t;
    mpz_t f;

    if (n % 2)
        abort();

    fmpq_init(b);
    mpfr_init2(t, mpfr_get_prec(x) + FLINT_BIT_COUNT(n) + 2);
    mpz_init(f);

    bernoulli_number(b, n);

    fmpq_get_mpfr(x, b, MPFR_RNDD);
    mpfr_const_pi(t, MPFR_RNDD);
    mpfr_mul_2exp(t, t, 1, MPFR_RNDD);
    mpfr_pow_ui(t, t, n, MPFR_RNDD);
    mpz_fac_ui(f, n);
    mpfr_div_z(t, t, f, MPFR_RNDD);
    mpfr_mul(x, x, t, MPFR_RNDD);
    mpfr_abs(x, x, MPFR_RNDD);
    mpfr_div_2exp(x, x, 1, MPFR_RNDD);

    mpfr_clear(t);
    fmpq_clear(b);
    mpz_clear(f);
}

void
arb_zeta_ui(arb_t x, ulong n)
{
    long prec = arb_prec(x);

    /* asymptotic case */
    if (n == 0 || n > 0.7 * prec)
    {
        arb_zeta_ui_mpfr(x, n);
    }
    else if (n == 3)
    {
        arb_const_zeta3_bsplit(x);
    }
    /* small even n */
    else if ((n % 2 == 0) && (n < 40 + 0.11*prec))
    {
        mpfr_t t;
        mpfr_init2(t, prec + 10);
        mpfr_zeta_ui_bernoulli(t, n);
        arb_set_mpfr(x, t, 100);  /* XXX */
        mpfr_clear(t);
    }
    /* small odd n, extremely high precision */
    /* FIXME: arb: n < prec * 0.0006 */
    else if (n < prec * 0.0006)
    {
        arb_zeta_ui_bsplit(x, n);
    }
    /* large n */
    else if (prec > 20 && n > 0.4 * pow(prec, 0.8))
    {
        arb_zeta_ui_euler_product(x, n);
    }
    /* fallback */
    else
    {
        arb_zeta_ui_mpfr(x, n);
    }
}
