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

/* tuning factor */
#define GAMMA_STIRLING_BETA 0.27

static __inline__ long
_fmpr_mag(const fmpr_t x)
{
    if (fmpr_is_special(x))
    {
        if (fmpr_is_zero(x))
            return -FMPR_PREC_EXACT;
        else
            return FMPR_PREC_EXACT;
    }
    else
    {
        return fmpz_bits(fmpr_manref(x)) + *fmpr_expref(x) - 1;
    }
}

#define PI 3.1415926535897932385

static long
choose_n(double log2z, double argz, int digamma, long prec)
{
    double argf, boundn;
    long n;

    argf = 1.0 / cos(0.5 * argz);
    argf = log(argf) * (1. / log(2));

    for (n = 1; ; n++)
    {
        if (digamma)
            boundn = bernoulli_bound_2exp_si(2*n) - (2*n)*log2z + (2*n+1)*argf;
        else
            boundn = bernoulli_bound_2exp_si(2*n) - (2*n-1)*log2z + (2*n)*argf;

        /* success */
        if (boundn <= -prec)
            return n;

        /* if the term magnitude does not decrease, r is too small */
        if (boundn > 1)
        {
            printf("exception: gamma_stirling_choose_param failed to converge\n");
            abort();
        }
    }
}

void
choose_small(int * reflect, long * r, long * n,
    double x, double y, int use_reflect, int digamma, long prec)
{
    double w, argz, log2z;
    long rr;

    /* use reflection formula if very negative */
    if (x < -5.0 && use_reflect)
    {
        *reflect = 1;
        x = 1.0 - x;
    }
    else
    {
        *reflect = 0;
    }

    /* argument reduction until |z| >= w */
    w = FLINT_MAX(1.0, GAMMA_STIRLING_BETA * prec);

    rr = 0;
    while (x < 1.0 || x*x + y*y < w*w)
    {
        x++;
        rr++;
    }

    log2z = 0.5 * log(x*x + y*y) * 1.44269504088896341;
    argz = atan2(y, x);

    *r = rr;
    *n = choose_n(log2z, argz, digamma, prec);
}

void
choose_large(int * reflect, long * r, long * n,
    const fmpr_t a, const fmpr_t b, int use_reflect, int digamma, long prec)
{
    if (use_reflect && fmpr_sgn(a) < 0)
        *reflect = 1;
    else
        *reflect = 0;

    *r = 0;

    /* so big that we will certainly have n = 0 */
    if (fmpr_cmpabs_2exp_si(a, LONG_MAX / 8) >= 0 ||
        fmpr_cmpabs_2exp_si(b, LONG_MAX / 8) >= 0)
    {
        *n = 0;
    }
    else
    {
        long ab, bb;
        double log2z, argz;

        ab = _fmpr_mag(a);
        bb = _fmpr_mag(b);

        log2z = FLINT_MAX(ab, bb);

        /* piecewise approximation of the argument */
        if (fmpr_is_zero(b))
        {
            if ((fmpr_sgn(a) < 0) && !(*reflect))
                argz = PI;
            else
                argz = 0.0;
        }
        else
        {
            if ((fmpr_sgn(a) < 0) && !(*reflect))
                if (fmpr_cmpabs(a, b) <= 0)
                    argz = PI * 0.75;
                else
                    argz = PI;
            else
                if (fmpr_cmpabs(a, b) <= 0)
                    argz = PI * 0.25;
                else
                    argz = PI * 0.5;
        }

        if (argz == PI)
            *n = 0;
        else
            *n = choose_n(log2z, argz, digamma, prec);
    }
}


void
gamma_stirling_choose_param_fmpcb(int * reflect, long * r, long * n,
    const fmpcb_t z, int use_reflect, int digamma, long prec)
{
    const fmpr_struct * a = fmprb_midref(fmpcb_realref(z));
    const fmpr_struct * b = fmprb_midref(fmpcb_imagref(z));

    if (fmpr_is_inf(a) || fmpr_is_nan(a) || fmpr_is_inf(b) || fmpr_is_nan(b))
    {
        *reflect = *r = *n = 0;
    }
    else if (fmpr_cmpabs_2exp_si(a, 40) > 0 || fmpr_cmpabs_2exp_si(b, 40) > 0)
    {
        choose_large(reflect, r, n, a, b, use_reflect, digamma, prec);
    }
    else
    {
        choose_small(reflect, r, n,
            fmpr_get_d(a, FMPR_RND_NEAR),
            fmpr_get_d(b, FMPR_RND_NEAR), use_reflect, digamma, prec);
    }
}

void
gamma_stirling_choose_param_fmprb(int * reflect, long * r, long * n,
    const fmprb_t x, int use_reflect, int digamma, long prec)
{
    const fmpr_struct * a = fmprb_midref(x);

    if (fmpr_is_inf(a) || fmpr_is_nan(a))
    {
        *reflect = *r = *n = 0;
    }
    else if (fmpr_cmpabs_2exp_si(a, 40) > 0)
    {
        fmpr_t b;
        fmpr_init(b);
        choose_large(reflect, r, n, a, b, use_reflect, digamma, prec);
        fmpr_clear(b);
    }
    else
    {
        choose_small(reflect, r, n,
            fmpr_get_d(a, FMPR_RND_NEAR), 0.0, use_reflect, digamma, prec);
    }
}

