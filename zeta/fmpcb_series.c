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

#include "zeta.h"

void
zeta_series(fmpcb_ptr z, const fmpcb_t s, const fmpcb_t a, int deflate, long d, long prec)
{
    ulong M, N;
    long i;
    fmpr_t bound;
    fmprb_ptr vb;

    if (d < 1)
        return;

    fmpr_init(bound);
    vb = _fmprb_vec_init(d);

    zeta_series_em_choose_param(bound, &N, &M, s, a, d, prec, FMPRB_RAD_PREC);
    zeta_series_em_vec_bound(vb, s, a, N, M, d, FMPRB_RAD_PREC);

    zeta_series_em_sum(z, s, a, deflate, N, M, d, prec);

    for (i = 0; i < d; i++)
    {
        fmprb_get_abs_ubound_fmpr(bound, vb + i, FMPRB_RAD_PREC);
        fmprb_add_error_fmpr(fmpcb_realref(z + i), bound);
        fmprb_add_error_fmpr(fmpcb_imagref(z + i), bound);
    }

    fmpr_clear(bound);
    _fmprb_vec_clear(vb, d);
}

