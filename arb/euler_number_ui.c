/*
    Copyright (C) 2016 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb.h"
#include "acb_dirichlet.h"

#if FLINT64
#define ARB_EULER_NUMBER_TAB_SIZE 25
#else
#define ARB_EULER_NUMBER_TAB_SIZE 15
#endif

static const ulong arb_euler_number_tab[] =
{
    1, 1, 5, 61, 1385, 50521, 2702765, 199360981,
#if FLINT64
    UWORD(19391512145), UWORD(2404879675441), UWORD(370371188237525),
    UWORD(69348874393137901), UWORD(15514534163557086905)
#endif
};

static double
arb_euler_number_mag(double n)
{
    double x;
    x = n + 2;
    x += ((n + 1) * log(n + 1) - n) * 1.44269504088897;  /* 1/log(2) */
    x -= 1.6514961294723*(n+1);  /* log2(pi) */
    return x;
}

static void
arb_euler_number_ui_beta(arb_t res, ulong n, slong prec)
{
    slong pi_prec;
    arb_t t;
    const signed char chi[4] = {0, 1, 0, -1};

    pi_prec = prec + 2 * FLINT_BIT_COUNT(n);
    arb_init(t);

    /* |E_n| = 2 n! beta(n+1) / (pi/2)^(n+1) */
    arb_const_pi(t, pi_prec);
    arb_mul_2exp_si(t, t, -1);
    arb_pow_ui(t, t, n + 1, pi_prec);

    _acb_dirichlet_euler_product_real_ui(res, n + 1, chi, 4, 1, prec);

    arb_mul(res, res, t, prec);
    arb_fac_ui(t, n, pi_prec);  /* todo: prec should be enough */
    arb_div(res, t, res, prec);
    arb_mul_2exp_si(res, res, 1);

    if (n % 4 == 2)
        arb_neg(res, res);

    arb_clear(t);
}

void
arb_euler_number_ui(arb_t res, ulong n, slong prec)
{
    double mag;

    if (n % 2)
    {
        arb_zero(res);
        return;
    }

    if (n < ARB_EULER_NUMBER_TAB_SIZE)
    {
        arb_set_ui(res, arb_euler_number_tab[n / 2]);
        if (n % 4 == 2)
            arb_neg(res, res);
        arb_set_round(res, res, prec);
        return;
    }

    mag = arb_euler_number_mag(n);

    if (prec > 0.9 * mag)
    {
        fmpz_t t;
        fmpz_init(t);
        arb_euler_number_ui_beta(res, n, mag + 5);
        if (arb_get_unique_fmpz(t, res))
            arb_set_round_fmpz(res, t, prec);
        fmpz_clear(t);
    }
    else
    {
        arb_euler_number_ui_beta(res, n, prec + 5);
        arb_set_round(res, res, prec);
    }
}

