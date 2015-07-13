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

    Copyright (C) 2015 Fredrik Johansson

******************************************************************************/

#include "acb_hypgeom.h"

/*

[S(k+1)] = [ R(k)  0   ] [S(k)]
[T(k+1)]   [ 1     1   ] [T(k)]

[S(k+1)] = [ P(k) / Q(k)  0   ] [S(k)]
[T(k+1)]   [ 1            1   ] [T(k)]


  1  [ P(k)         ]
---- [              ]
Q(k) [ Q(k)   Q(k)  ]

[[A2 0] [B2 C2]] . [[A1 0] [B1 C1]] = [[A1 A2 0] [A1 B2 + B1 C2   C1 C2]

A1 B2 + B1 B2 = B2 (A1 + B1) -- use to save time?

*/

static void
bsplit(acb_t A1, acb_t B1, acb_t C1,
        acb_srcptr a, long p,
        acb_srcptr b, long q,
        const acb_t z,
        long aa,
        long bb,
        long prec,
        int invz)
{
    if (bb - aa == 1)
    {
        long i;

        if (p == 0)
        {
            if (invz)
                acb_one(A1);
            else
                acb_set(A1, z);
        }
        else
        {
            acb_add_ui(A1, a, aa, prec);

            for (i = 1; i < p; i++)
            {
                acb_add_ui(B1, a + i, aa, prec);
                acb_mul(A1, A1, B1, prec);
            }

            if (!invz)
                acb_mul(A1, A1, z, prec);
        }

        if (q == 0)
        {
            if (invz)
                acb_set(C1, z);
            else
                acb_one(C1);
        }
        else
        {
            acb_add_ui(C1, b, aa, prec);

            for (i = 1; i < q; i++)
            {
                acb_add_ui(B1, b + i, aa, prec);
                acb_mul(C1, C1, B1, prec);
            }

            if (invz)
                acb_mul(C1, C1, z, prec);
        }

        /* acb_set(B1, C1);   but we skip this */
    }
    else
    {
        long m;

        acb_t A2, B2, C2;

        acb_init(A2);
        acb_init(B2);
        acb_init(C2);

        m = aa + (bb - aa) / 2;

        bsplit(A1, B1, C1, a, p, b, q, z, aa, m, prec, invz);
        bsplit(A2, B2, C2, a, p, b, q, z, m, bb, prec, invz);

        if (bb - m == 1)  /* B2 = C2 */
        {
            if (m - aa == 1)
                acb_add(B2, A1, C1, prec);
            else
                acb_add(B2, A1, B1, prec);

            acb_mul(B1, B2, C2, prec);
        }
        else
        {
            if (m - aa == 1)
                acb_mul(B1, C1, C2, prec);
            else
                acb_mul(B1, B1, C2, prec);

            acb_addmul(B1, A1, B2, prec);
        }

        acb_mul(A1, A1, A2, prec);
        acb_mul(C1, C1, C2, prec);

        acb_clear(A2);
        acb_clear(B2);
        acb_clear(C2);
    }
}

void
acb_hypgeom_pfq_sum_bs(acb_t s, acb_t t,
    acb_srcptr a, long p, acb_srcptr b, long q, const acb_t z, long n, long prec)
{
    acb_t u, v, w;

    if (n < 4)
    {
        acb_hypgeom_pfq_sum_forward(s, t, a, p, b, q, z, n, prec);
        return;
    }

    acb_init(u);
    acb_init(v);
    acb_init(w);

    bsplit(u, v, w, a, p, b, q, z, 0, n, prec, 0);

    acb_div(t, u, w, prec);
    acb_div(s, v, w, prec);

    acb_clear(u);
    acb_clear(v);
    acb_clear(w);
}

void
acb_hypgeom_pfq_sum_bs_invz(acb_t s, acb_t t,
    acb_srcptr a, long p, acb_srcptr b, long q, const acb_t z, long n, long prec)
{
    acb_t u, v, w;

    if (n < 4)
    {
        acb_init(u);
        acb_inv(u, z, prec);
        acb_hypgeom_pfq_sum_forward(s, t, a, p, b, q, u, n, prec);
        acb_clear(u);
        return;
    }

    acb_init(u);
    acb_init(v);
    acb_init(w);

    bsplit(u, v, w, a, p, b, q, z, 0, n, prec, 1);

    acb_div(t, u, w, prec);
    acb_div(s, v, w, prec);

    acb_clear(u);
    acb_clear(v);
    acb_clear(w);
}

