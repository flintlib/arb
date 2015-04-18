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

#include "acb_poly.h"

static ulong choose_M(ulong N, long target)
{
    return FLINT_MIN(N, target + N / 100);
}

void
_acb_poly_zeta_em_choose_param(arf_t bound, ulong * N, ulong * M, const acb_t s, const acb_t a, long d, long target, long prec)
{
    ulong A, B, C, limit;
    arf_t Abound, Bbound, Cbound, tol;

    arf_init(Abound);
    arf_init(Bbound);
    arf_init(Cbound);
    arf_init(tol);

    arf_set_si_2exp_si(tol, 1, -target);

    A = 1;
    B = 2;

    /* Hack: allow evaluation very high up in the critical strip... */
    if (arf_cmpabs_2exp_si(arb_midref(acb_imagref(s)), 10) > 0)
        limit = ULONG_MAX / 4;
    else
        limit = 100 * target;  /* but normally, fail more quickly */

    _acb_poly_zeta_em_bound1(Bbound, s, a, B, choose_M(B, target), d, prec);

    if (arf_cmp(Bbound, tol) > 0)
    {
        while (arf_cmp(Bbound, tol) > 0 && B <= limit)
        {
            arf_set(Abound, Bbound);
            A *= 2;
            B *= 2;

            if (B == 0) abort();

            _acb_poly_zeta_em_bound1(Bbound, s, a, B, choose_M(B, target), d, prec);
        }

        /* bisect (-A, B] */
        while (B > A + 4)
        {
            C = A + (B - A) / 2;

            _acb_poly_zeta_em_bound1(Cbound, s, a, C, choose_M(C, target), d, prec);

            if (arf_cmp(Cbound, tol) < 0)
            {
                B = C;
                arf_set(Bbound, Cbound);
            }
            else
            {
                A = C;
                arf_set(Abound, Cbound);
            }
        }
    }

    arf_set(bound, Bbound);
    *N = B;
    *M = choose_M(B, target);

    arf_clear(Abound);
    arf_clear(Bbound);
    arf_clear(Cbound);
    arf_clear(tol);
}

