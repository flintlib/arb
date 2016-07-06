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
acb_dirichlet_number_primitive(const acb_dirichlet_group_t G)
{
    if (G->q % 4 == 2)
        return 0;
    else
    {
        slong k;
        ulong n = 1;

        /* no overflow since result < G->q */

        for (k = (G->neven == 2); k < G->num; k++)
        {
            ulong p = G->P[k].p, e = G->P[k].e;
            if (e == 1)
                n *= p - 2;
            else
                n *= (p * (p - 2) + 1) * n_pow(p, e - 2);
        }

        return n;
    }
}


