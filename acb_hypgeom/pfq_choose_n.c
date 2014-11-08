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

long
acb_hypgeom_pfq_choose_n(acb_srcptr a, long p,
                         acb_srcptr b, long q, const acb_t z, long prec)
{
    long k, n, minimum_n, maximum_n;
    mag_t zmag;

    /* todo: first test for z = 0 and nonpositive integers */

    double t, u;
    double log2_z;
    double log2_term;
    double log2_factor;
    double log2_term_max;

    double * are;
    double * aim;
    double * bre;
    double * bim;

    mag_init(zmag);

    are = flint_malloc(sizeof(double) * 2 * (p + q));
    aim = are + p;
    bre = aim + p;
    bim = bre + q;

    acb_get_mag(zmag, z);
    log2_z = mag_get_log2_d_approx(zmag);

    for (k = 0; k < p; k++)
    {
        are[k] = arf_get_d(arb_midref(acb_realref(a + k)), ARF_RND_DOWN);
        aim[k] = arf_get_d(arb_midref(acb_imagref(a + k)), ARF_RND_DOWN);
    }

    minimum_n = 1;
    maximum_n = 50 + 2 * prec;

    for (k = 0; k < q; k++)
    {
        bre[k] = arf_get_d(arb_midref(acb_realref(b + k)), ARF_RND_DOWN);
        bim[k] = arf_get_d(arb_midref(acb_imagref(b + k)), ARF_RND_DOWN);

        if (bre[k] <= 0.25)
        {
            minimum_n = FLINT_MAX(minimum_n, 2 - bre[k]);
        }
    }

    n = 1;
    log2_term = 0.0;
    log2_term_max = log2_term;

    while (n <= maximum_n && minimum_n < maximum_n)
    {
        if (log2_term < log2_term_max - prec - 4 && n >= minimum_n)
            break;

        t = 1.0;

        for (k = 0; k < FLINT_MAX(p, q); k++)
        {
            if (k < p)
            {
                u = (are[k] + n) * (are[k] + n) + (aim[k] * aim[k]);
                u = FLINT_ABS(u);

                if (u < 1e-8 || u > 1e100 || t > 1e100)
                    goto somethingstrange;

                t *= u;
            }

            if (k < q)
            {
                u = (bre[k] + n) * (bre[k] + n) + (bim[k] * bim[k]);
                u = FLINT_ABS(u);

                if (u < 1e-8 || u > 1e100 || t > 1e100)
                    goto somethingstrange;

                t /= u;
            }
        }

        log2_factor = 0.5 * log(t) * 1.4426950408889634074 + log2_z;

        /* For asymptotic series, require rapid decay */
        if (p > q && n >= minimum_n && log2_factor > -0.2)
            break;

        log2_term += log2_factor;
        log2_term_max = FLINT_MAX(log2_term_max, log2_term);
        n++;
    }

somethingstrange:

    flint_free(are);
    mag_clear(zmag);

    return n;
}

