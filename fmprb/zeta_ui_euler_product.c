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

#include <math.h>
#include "fmprb.h"

/*
Let P(a,b) = prod_{a <= p <= b} (1 - p^(-s)).
Then 1/zeta(s) = P(a,M) * P(M+1,inf).

According to the analysis in S. Fillebrown,
"Faster Computation of Bernoulli Numbers", Journal of Algorithms 13,
431-445 (1992), it holds for all s >= 6 and M >= 1 that
(1/P(M+1,inf) - 1) <= 2 * M^(1-s) / (s/2 - 1).

Writing 1/zeta(s) = P(a,M) * (1 - eps) and solving for eps gives

1/(1-eps) <= 1 + 2 * M^(1-s) / (s/2 - 1), so we have
eps <= 2 * M^(1-s) / (s/2 - 1).

Since 0 < P(a,M) <= 1, this bounds the absolute error of 1/zeta(s).
*/

void
fmprb_zeta_inv_ui_euler_product(fmprb_t z, ulong s, long prec)
{
    long wp, powprec;
    fmprb_t t;
    mp_limb_t p;

    if (s < 6)
    {
        printf("too small s!\n");
        abort();
    }

    /* heuristic */
    wp = prec + FLINT_BIT_COUNT(prec) + (prec/s) + 4;

    fmprb_init(t);

    /* z = 1 */
    fmprb_set_ui(z, 1UL);

    /* z = 1 - 2^(-s) */
    {
        fmprb_t w;
        fmprb_init(w);
        fmprb_set_ui(w, 1);
        fmpz_set_si(fmpr_expref(fmprb_midref(w)), -s);
        fmprb_sub(z, z, w, wp);
        fmprb_clear(w);
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

        fmprb_ui_pow_ui(t, p, s, powprec);
        fmprb_div(t, z, t, powprec);
        fmprb_sub(z, z, t, wp);

        p = n_nextprime(p, 0);
    }

    /* Truncation error based on the termination test */
    fmprb_add_error_2exp_si(z, -(prec+1));

    fmprb_clear(t);

    /* TODO: change precision to prec here */
}

void
fmprb_zeta_ui_euler_product(fmprb_t z, ulong s, long prec)
{
    fmprb_t one;
    fmprb_init(one);
    fmprb_set_ui(one, 1);
    fmprb_zeta_inv_ui_euler_product(z, s, prec);
    fmprb_div(z, one, z, prec);
    fmprb_clear(one);
}
