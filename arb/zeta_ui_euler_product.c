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
#include "arb.h"

static __inline__ void
arb_mul_ui(arb_t y, const arb_t x, ulong c)
{
    fmpz_mul_ui(arb_midref(y), arb_midref(x), c);
    fmpz_mul_ui(arb_radref(y), arb_radref(x), c);
    fmpz_set(arb_expref(y), arb_expref(x));
}

static __inline__ void
arb_ui_pow_ui(arb_t y, ulong b, ulong e)
{
    long i;

    if (e <= 1)
    {
        arb_set_ui(y, e == 0 ? 1 : b);
        return;
    }

    arb_set_ui(y, b);
    for (i = FLINT_BIT_COUNT(e) - 2; i >= 0; i--)
    {
        arb_mul(y, y, y);
        if (e & (1UL<<i))
            arb_mul_ui(y, y, b);
    }
}

void
arb_zeta_inv_ui_euler_product(arb_t z, ulong s)
{
    long prec, wp, powprec;
    arb_t t;
    mp_limb_t p;

    if (s < 6)
    {
        printf("too small s!\n");
        abort();
    }

    prec = arb_prec(z);
    wp = prec + FLINT_BIT_COUNT(prec) + (prec/s) + 4;

    arb_init(t, wp);
    z->prec = wp;

    /* z = 1 */
    arb_set_ui(z, 1UL);
    fmpz_mul_2exp(arb_midref(z), arb_midref(z), wp);
    fmpz_set_si(arb_expref(z), -wp);

    /* z = 1 - 2^(-s) */
    if (wp >= s)
    {
        fmpz_t r;
        fmpz_init(r);
        fmpz_set_ui(r, 1UL);
        fmpz_mul_2exp(r, r, wp - s);
        fmpz_sub(arb_midref(z), arb_midref(z), r);
        fmpz_clear(r);
    }
    else
    {
        fmpz_set_ui(arb_radref(z), 1);
    }

    p = 3UL;

    while (1)
    {
        /* approximate magnitude of p^s */
        double powmag = s * log(p) * 1.4426950408889634;
        powprec = FLINT_MAX(wp - powmag, 8);

        /* see error analysis */
        if ((powmag >= prec) &&
            -((s-1)*log(p-1)) - log(s/2-1) + 1 <= -(prec+1) * 0.69314718055995)
                break;

        arb_prec(t) = powprec;
        arb_ui_pow_ui(t, p, s);
        arb_div(t, z, t);
        arb_sub(z, z, t);

        p = n_nextprime(p, 0);
    }

    /* Truncation error based on the termination test */
    arb_add_error_2exp(z, -(prec+1));

    arb_clear(t);
    z->prec = prec;
    _arb_normalise(z);
}

void
arb_zeta_ui_euler_product(arb_t z, ulong s)
{
    arb_t one;
    arb_init(one, 1);
    arb_set_ui(one, 1);
    arb_zeta_inv_ui_euler_product(z, s);
    arb_div(z, one, z);
    arb_clear(one);
}
