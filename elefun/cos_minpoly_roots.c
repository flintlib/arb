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

#include "elefun.h"
#include "ulong_extras.h"

/* computes all the d roots alpha_j of Phi_n(x) */
void
_elefun_cos_minpoly_roots(fmprb_ptr alpha, long d, ulong n, long prec)
{
    fmprb_t t, u, v, s1, s2, c1, c2;
    long i, j;

    fmprb_init(t);
    fmprb_init(u);
    fmprb_init(v);
    fmprb_init(s1);
    fmprb_init(s2);
    fmprb_init(c1);
    fmprb_init(c2);

    fmprb_const_pi(s1, prec);
    fmprb_mul_ui(s1, s1, 2, prec);
    fmprb_div_ui(s1, s1, n, prec);
    fmprb_sin_cos(s1, c1, s1, prec);
    fmprb_one(c2);
    fmprb_zero(s2);

    for (i = j = 0; j < d; i++)
    {
        if (n_gcd(n, i) == 1)
        {
            fmprb_set(alpha + j, c2);
            j++;
        }

        fmprb_mul(t, c1, c2, prec);
        fmprb_mul(u, s1, s2, prec);
        fmprb_sub(v, t, u, prec);

        fmprb_mul(t, c1, s2, prec);
        fmprb_mul(u, s1, c2, prec);
        fmprb_add(s2, t, u, prec);
        fmprb_set(c2, v);
    }

    fmprb_clear(t);
    fmprb_clear(u);
    fmprb_clear(v);
    fmprb_clear(s1);
    fmprb_clear(s2);
    fmprb_clear(c1);
    fmprb_clear(c2);
}

