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

    Copyright (C) 2015 Jonathan Bober
    Copyright (C) 2016 Fredrik Johansson

******************************************************************************/

#include "acb_dirichlet.h"

void
acb_dirichlet_acb_pairing(acb_t res, const acb_dirichlet_group_t G, ulong m, ulong n, slong prec)
{
    ulong expo;
    expo = acb_dirichlet_pairing(G, m, n);
    if (expo == ACB_DIRICHLET_CHI_NULL)
        acb_zero(res);
    else
    {
        fmpq_t t;
        fmpq_init(t);
        fmpq_set_si(t, 2 * expo, G->expo);
        arb_sin_cos_pi_fmpq(acb_imagref(res), acb_realref(res), t, prec);
        fmpq_clear(t);
    }
}
