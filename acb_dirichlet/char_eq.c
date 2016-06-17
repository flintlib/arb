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

int
acb_dirichlet_char_eq(const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi1, const acb_dirichlet_char_t chi2)
{
    acb_dirichlet_conrey_t x, y;
    
    if (chi1->q != chi2->q)
        return 0;

    if (chi1->order != chi2->order)
        return 0;

    if (chi1->conductor != chi2->conductor)
        return 0;

    if (!acb_dirichlet_conrey_eq(G, chi1->x, chi2->x))
        return 0;

    x->n = y->n = 1;
    x->log = chi1->expo;
    y->log = chi2->expo;
    if (!acb_dirichlet_conrey_eq(G, x, y))
        return 0;

    return 1;
}
