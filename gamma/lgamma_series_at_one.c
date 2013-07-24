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
#include "zeta.h"

void
gamma_lgamma_series_at_one(fmprb_ptr u, long len, long prec)
{
    long i;

    if (len > 0) fmprb_zero(u);
    if (len > 1) fmprb_const_euler(u + 1, prec);
    if (len > 2) zeta_ui_vec(u + 2, 2, len - 2, prec);

    for (i = 2; i < len; i++)
        fmprb_div_ui(u + i, u + i, i, prec);

    for (i = 1; i < len; i += 2)
        fmprb_neg(u + i, u + i);
}

