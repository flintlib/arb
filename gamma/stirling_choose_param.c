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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "gamma.h"
#include "bernoulli.h"

double
gamma_stirling_arg_factor(double x, double y)
{
    double t;

    if (y == 0)
    {
        if (x <= 0)
            return 1e30;
        else
            return 0.0;
    }
    else
    {
        t = 1.0 / cos(0.5 * atan2(y, x));
        return 2.0 * log(t) * (1. / log(2));
    }
}

void
gamma_stirling_choose_param(int * reflect, long * r, long * n,
    double x, double y, double beta, int allow_reflection, long prec)
{
    double w, logz, argz, boundn, prevboundn;
    int i;
    long rr, nn;

    /* use reflection formula if very negative */
    if (x < -5.0 && allow_reflection)
    {
        *reflect = 1;
        x = 1.0 - x;
    }
    else
    {
        *reflect = 0;
    }

    /* argument reduction until |z| >= w */
    w = FLINT_MAX(1.0, beta * prec);
    rr = 0;
    while (x < 1.0 || x*x + y*y < w*w)
    {
        x++;
        rr++;
    }

    /* this should succeed on the first try if the selection strategy is sound
       (i.e. if GAMMA_STIRLING_BETA is large enough) */
    for (i = 0; i < 1; i++)
    {
        /* log2(z) */
        logz = 0.5 * log(x*x + y*y) * 1.44269504088896341;

        /* log2(1/cos(0.5*arg(z))^2) */
        argz = gamma_stirling_arg_factor(x, y);

        nn = 1;
        prevboundn = 1e300;

        while (1)
        {
            boundn = bernoulli_bound_2exp_si(2 * nn) - (2 * nn-1) * logz + argz * nn;

            /* success */
            if (boundn <= -prec)
            {
                *r = rr;
                *n = nn;
                return;
            }

            /* if the term magnitude does not decrease, try larger r */
            if (boundn > 1)
            {
                int r2 = 4 + rr * 0.5;
                x += r2;
                rr += r2;
                break;
            }

            prevboundn = boundn;
            nn++;
        }
    }

    printf("exception: gamma_stirling_choose_param failed to converge\n");
    abort();
}

