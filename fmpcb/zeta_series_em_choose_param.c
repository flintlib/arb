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

#include "fmpcb.h"

static ulong choose_M(ulong N, long target)
{
    return FLINT_MIN(N, target + N / 100);
}

void
fmpcb_zeta_series_em_choose_param(fmpr_t bound, ulong * N, ulong * M, const fmpcb_t s, const fmpcb_t a, long d, long target, long prec)
{
    ulong A, B, C;
    fmpr_t Abound, Bbound, Cbound, tol;

    fmpr_init(Abound);
    fmpr_init(Bbound);
    fmpr_init(Cbound);
    fmpr_init(tol);

    fmpr_set_si_2exp_si(tol, 1, -target);

    A = 1;
    B = 2;

    fmpcb_zeta_series_em_bound(Bbound, s, a, B, choose_M(B, target), d, prec);

    if (fmpr_cmp(Bbound, tol) > 0)
    {
        while (fmpr_cmp(Bbound, tol) > 0)
        {
            fmpr_set(Abound, Bbound);
            A *= 2;
            B *= 2;

            if (B == 0) abort();

            fmpcb_zeta_series_em_bound(Bbound, s, a, B, choose_M(B, target), d, prec);
        }

        /* bisect (-A, B] */
        while (B > A + 4)
        {
            C = A + (B - A) / 2;

            fmpcb_zeta_series_em_bound(Cbound, s, a, C, choose_M(C, target), d, prec);

            if (fmpr_cmp(Cbound, tol) < 0)
            {
                B = C;
                fmpr_set(Bbound, Cbound);
            }
            else
            {
                A = C;
                fmpr_set(Abound, Cbound);
            }
        }
    }

    fmpr_set(bound, Bbound);
    *N = B;
    *M = choose_M(B, target);

    fmpr_clear(Abound);
    fmpr_clear(Bbound);
    fmpr_clear(Cbound);
    fmpr_clear(tol);
}
