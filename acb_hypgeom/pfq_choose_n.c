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

#include "acb_hypgeom.h"

double mag_get_log2_d_approx(const mag_t x);

int
acb_hypgeom_pfq_choose_n_double(long * nn,
    const double * are, const double * aim, long p,
    const double * bre, const double * bim, long q,
    double log2_z,
    long n_skip, long n_min, long n_max, long prec)
{
    double increase, term, term_max, accuracy, accuracy_best, t, u;
    double required_decrease;
    long k, n, n_best;
    int success;

    if (p < q)
        required_decrease = 0.01;
    else if (p == q)
        required_decrease = 0.0001;
    else
        required_decrease = 0.01;

    term = term_max = accuracy = accuracy_best = 0.0;
    success = 0;

    for (n = n_best = n_skip; n < n_max; n++)
    {
        t = 1.0;

        for (k = 0; k < FLINT_MAX(p, q); k++)
        {
            if (k < p)
            {
                u = (are[k]+n-1)*(are[k]+n-1) + (aim[k]*aim[k]);
                t *= FLINT_ABS(u);
            }

            if (k < q)
            {
                u = (bre[k]+n-1)*(bre[k]+n-1) + (bim[k]*bim[k]);
                u = FLINT_ABS(u);

                if (u > 1e-100)
                    t /= u;
            }
        }

        increase = 0.5 * log(t) * 1.4426950408889634074 + log2_z;

        term += increase;
        term_max = FLINT_MAX(term_max, term);
        accuracy = term_max - term;

        if (accuracy > accuracy_best && n >= n_min && increase < -required_decrease)
        {
            n_best = n;
            accuracy_best = accuracy;
        }

        if (accuracy_best > prec + 4)
        {
            success = 1;
            break;
        }
    }

    *nn = n_best;
    return success;
}

long
acb_hypgeom_pfq_choose_n(acb_srcptr a, long p,
                         acb_srcptr b, long q, const acb_t z, long prec)
{
    long n_skip, n_min, n_max, n_terminating, nint;
    long k, n;
    double log2_z;
    double * are;
    double * aim;
    double * bre;
    double * bim;
    mag_t zmag;
    int success;

    if (acb_is_zero(z) || !acb_is_finite(z))
        return 1;

    mag_init(zmag);

    are = flint_malloc(sizeof(double) * 2 * (p + q));
    aim = are + p;
    bre = aim + p;
    bim = bre + q;

    acb_get_mag(zmag, z);
    log2_z = mag_get_log2_d_approx(zmag);

    n_skip = 1;
    n_min = 1;
    n_max = FLINT_MIN(LONG_MAX / 2, 50 + 10.0 * prec);
    n_terminating = LONG_MAX;

    for (k = 0; k < p; k++)
    {
        are[k] = arf_get_d(arb_midref(acb_realref(a + k)), ARF_RND_DOWN);
        aim[k] = arf_get_d(arb_midref(acb_imagref(a + k)), ARF_RND_DOWN);

        /* If the series is terminating, stop at this n. */
        if (acb_is_int(a + k) && are[k] <= 0.0)
        {
            n_terminating = FLINT_MIN(n_terminating, (long) (-are[k] + 1));
            n_terminating = FLINT_MAX(n_terminating, 1);
        }
        else if (are[k] <= 0.01 && FLINT_ABS(aim[k]) < 0.01)
        {
            /* If very close to an integer, we don't want to deal with it
               using doubles, so fast forward (todo: could work with the
               log2 of the difference instead when this happens). */
            nint = floor(are[k] + 0.5);
            if (FLINT_ABS(nint - are[k]) < 0.01)
                n_skip = FLINT_MAX(n_skip, 2 - nint);
        }
    }

    n_max = FLINT_MIN(n_max, n_terminating);

    for (k = 0; k < q; k++)
    {
        bre[k] = arf_get_d(arb_midref(acb_realref(b + k)), ARF_RND_DOWN);
        bim[k] = arf_get_d(arb_midref(acb_imagref(b + k)), ARF_RND_DOWN);

        if (bre[k] <= 0.25)
        {
            n_min = FLINT_MAX(n_min, 2 - bre[k]);

            /* Also avoid near-integers here (can even allow exact
               integers when computing regularized hypergeometric functions). */
            if (bre[k] <= 0.01 && FLINT_ABS(bim[k]) < 0.01)
            {
                nint = floor(bre[k] + 0.5);
                if (FLINT_ABS(nint - bre[k]) < 0.01)
                    n_skip = FLINT_MAX(n_skip, 2 - nint);
            }
        }
    }

    success = acb_hypgeom_pfq_choose_n_double(&n, are, aim, p, bre, bim, q,
        log2_z, n_skip, n_min, n_max, prec);

    if (!success)
    {
        if (n_terminating <= n_max)
        {
            n = n_terminating;
        }
        else
        {
            n = FLINT_MAX(n_min, n);
            n = FLINT_MIN(n_max, n);
        }
    }

    flint_free(are);
    mag_clear(zmag);

    return n;
}

long
acb_hypgeom_pfq_series_choose_n(const acb_poly_struct * a, long p,
                                const acb_poly_struct * b, long q,
                                const acb_poly_t z, long len, long prec)
{
    long n_skip, n_min, n_max, n_terminating, nint;
    long k, n;
    double log2_z;
    double * are;
    double * aim;
    double * bre;
    double * bim;
    acb_t t;
    mag_t zmag;
    int success;

    if (z->length == 0) /* todo: check finite */
        return 1;

    mag_init(zmag);
    acb_init(t);

    are = flint_malloc(sizeof(double) * 2 * (p + q));
    aim = are + p;
    bre = aim + p;
    bim = bre + q;

    acb_get_mag(zmag, z->coeffs);
    log2_z = mag_get_log2_d_approx(zmag);

    n_skip = 1;
    n_min = 1;
    n_max = FLINT_MIN(LONG_MAX / 2, 50 + 10.0 * prec);
    n_terminating = LONG_MAX;

    /* e.g. for exp(x + O(x^100)), make sure that we actually
       run the recurrence up to order 100.
       todo: compute correctly based on degrees of inputs */
    n_min = FLINT_MAX(n_min, len);
    n_max = FLINT_MAX(n_max, n_min);

    for (k = 0; k < p; k++)
    {
        acb_poly_get_coeff_acb(t, a + k, 0);

        are[k] = arf_get_d(arb_midref(acb_realref(t)), ARF_RND_DOWN);
        aim[k] = arf_get_d(arb_midref(acb_imagref(t)), ARF_RND_DOWN);

        /* If the series is terminating, stop at this n. */
        if (acb_is_int(t) && are[k] <= 0.0 && acb_poly_length(a + k) <= 1)
        {
            n_terminating = FLINT_MIN(n_terminating, (long) (-are[k] + 1));
            n_terminating = FLINT_MAX(n_terminating, 1);
        }
        else if (are[k] <= 0.01 && FLINT_ABS(aim[k]) < 0.01)
        {
            /* If very close to an integer, we don't want to deal with it
               using doubles, so fast forward (todo: could work with the
               log2 of the difference instead when this happens). */
            nint = floor(are[k] + 0.5);
            if (FLINT_ABS(nint - are[k]) < 0.01)
                n_skip = FLINT_MAX(n_skip, 2 - nint);
        }
    }

    n_max = FLINT_MIN(n_max, n_terminating);

    for (k = 0; k < q; k++)
    {
        acb_poly_get_coeff_acb(t, b + k, 0);

        bre[k] = arf_get_d(arb_midref(acb_realref(t)), ARF_RND_DOWN);
        bim[k] = arf_get_d(arb_midref(acb_imagref(t)), ARF_RND_DOWN);

        if (bre[k] <= 0.25)
        {
            n_min = FLINT_MAX(n_min, 2 - bre[k]);

            /* Also avoid near-integers here (can even allow exact
               integers when computing regularized hypergeometric functions). */
            if (bre[k] <= 0.01 && FLINT_ABS(bim[k]) < 0.01)
            {
                nint = floor(bre[k] + 0.5);
                if (FLINT_ABS(nint - bre[k]) < 0.01)
                    n_skip = FLINT_MAX(n_skip, 2 - nint);
            }
        }
    }

    success = acb_hypgeom_pfq_choose_n_double(&n, are, aim, p, bre, bim, q,
        log2_z, n_skip, n_min, n_max, prec);

    if (!success)
    {
        if (n_terminating <= n_max)
        {
            n = n_terminating;
        }
        else
        {
            n = FLINT_MAX(n_min, n);
            n = FLINT_MIN(n_max, n);
        }
    }

    flint_free(are);
    acb_clear(t);
    mag_clear(zmag);

    return n;
}

