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
#include "zeta.h"

/*
Let P(a,b) = prod_{a <= p <= b} (1 - p^(-s)).
Then 1/zeta(s) = P(a,M) * P(M+1,inf).

According to the analysis in S. Fillebrown,
"Faster Computation of Bernoulli Numbers", Journal of Algorithms 13,
431-445 (1992), it holds for all s >= 6 and M >= 1 that
(1/P(M+1,inf) - 1) <= 2 * M^(1-s) / (s/2 - 1).

Writing 1/zeta(s) = P(a,M) * (1 - eps) and solving for eps gives

1/(1-eps) <= 1 + 2 * M^(1-s) / (s/2 - 1), so we have
eps <= 2 * M^(1-s) / (s/2 - 1) = 4 * M^(1-s) / (s-2).

Since 0 < P(a,M) <= 1, this bounds the absolute error of 1/zeta(s).
*/

static void
add_error(fmprb_t z, ulong M, ulong s)
{
    fmpr_t t;
    fmpr_init(t);
    fmpr_set_ui(t, M);
    fmpr_pow_sloppy_ui(t, t, s - 1, FMPRB_RAD_PREC, FMPR_RND_DOWN);
    fmpr_mul_ui(t, t, s - 2, FMPRB_RAD_PREC, FMPR_RND_DOWN);
    fmpr_ui_div(t, 4, t, FMPRB_RAD_PREC, FMPR_RND_UP);
    fmprb_add_error_fmpr(z, t);
    fmpr_clear(t);
}

void
fmprb_zeta_inv_ui_euler_product(fmprb_t z, ulong s, long prec)
{
    long wp, powprec;
    double powmag;
    fmprb_t t;
    ulong M;
    mp_limb_t p;

    if (s < 6)
    {
        printf("too small s!\n");
        abort();
    }

    /* heuristic */
    wp = prec + FLINT_BIT_COUNT(prec) + (prec/s) + 4;

    fmprb_init(t);

    /* z = 1 - 2^(-s) */
    fmprb_set_ui(z, 1UL);
    fmpr_set_ui_2exp_si(fmprb_midref(t), 1, -s);
    fmprb_sub(z, z, t, wp);

    M = 2;
    p = 3UL;
    while (1)
    {
        /* approximate magnitude of p^s */
        powmag = s * log(p) * 1.4426950408889634;
        powprec = FLINT_MAX(wp - powmag, 8);

        /* see error analysis */
        if ((powmag >= prec) &&
            ((1.-s)*log(M-1.)) - log(s-2.) + 2 <= -(prec+1) * 0.69314718055995)
                break;

        M = p;
        fmprb_ui_pow_ui(t, p, s, powprec);
        fmprb_div(t, z, t, powprec);
        fmprb_sub(z, z, t, wp);

        p = n_nextprime(p, 0);
    }

    add_error(z, M, s);
    fmprb_clear(t);
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

