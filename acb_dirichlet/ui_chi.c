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

    Copyright (C) 2016 Pascal Molin

******************************************************************************/

#include "acb_dirichlet.h"

ulong
acb_dirichlet_ui_chi(const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, ulong n)
{
    if (n_gcd(G->q, n) > 1)
    {
        return ACB_DIRICHLET_CHI_NULL;
    }
    else
    {
        ulong v = 0, k;
        acb_dirichlet_conrey_t x;
        acb_dirichlet_conrey_init(x, G);
        acb_dirichlet_conrey_log(x, G, n);
        for (k = 0; k < G->num; k++)
            v = (v + chi->expo[k] * x->log[k]) % chi->order;
        acb_dirichlet_conrey_clear(x);
        return v;
    }
}
