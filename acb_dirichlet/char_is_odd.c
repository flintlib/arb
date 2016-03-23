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
acb_dirichlet_char_is_odd(const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi)
{
    slong k, odd = 0;
    for (k = 0; k < G->num; k++)
    {
        if (k == 1 && G->neven == 2)
            continue;
        if (chi->expo[k] % 2 == 1)
            odd = 1 - odd;
    }
    return odd;
}
