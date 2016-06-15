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
acb_dirichlet_conductor_ui(const acb_dirichlet_group_t G, ulong a)
{
    slong k;
    ulong ap, cond;
    nmod_t pe;

    cond = 1;

    for (k = (G->neven == 2); k < G->num; k++)
    {
        ulong p, e;
        p = G->primes[k];
        e = G->exponents[k];
        nmod_init(&pe, G->primepowers[k]);
        ap = a % pe.n;
        if (ap == 1)
            continue;
        if (p == 2)
        {
            cond = 4;
            if (a % 4 == 3)
                ap = pe.n - ap;
        }
        else
        {
            cond *= p;
            ap = nmod_pow_ui(ap, p - 1, pe);
        }

        while (ap != 1)
        {
            cond *= p;
            ap = nmod_pow_ui(ap, p, pe);
        }

    }

    return cond;
}
