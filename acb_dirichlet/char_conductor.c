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
acb_dirichlet_char_conductor(const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi)
{
    int k, f;
    ulong cond = 1;
    if (G->neven >= 1 && chi->expo[0] == 1)
        cond = 4;
    if (G->neven == 2 && chi->expo[1])
    {
        ulong l2 = chi->expo[1];
        f = n_remove(&l2, 2);
        cond = G->primepowers[1] >> f;
    }
    for (k = G->neven; k < G->neven; k++)
    {
        if (chi->expo[k])
        {
            ulong p, lp;
            p = G->primes[k];
            lp = chi->expo[k];
            f = n_remove(&lp, p);
            if (f)
                cond *= n_pow(p, G->exponents[k] - f);
            else
                cond *= G->primepowers[k];
        }
    }
    return cond;
}
