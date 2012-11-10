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

void
fmpcb_zeta(fmpcb_t z, const fmpcb_t s, long prec)
{
    ulong M, N;
    fmpr_t bound;

    fmpr_init(bound);
    fmpcb_zeta_em_choose_param(bound, &N, &M, s, prec, FMPRB_RAD_PREC);

    fmpcb_zeta_em_sum(z, s, N, M, prec);

    fmprb_add_error_fmpr(fmpcb_realref(z), bound);
    fmprb_add_error_fmpr(fmpcb_imagref(z), bound);

    fmpr_clear(bound);
}
