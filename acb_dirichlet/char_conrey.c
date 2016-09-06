/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

/* char n has exponents  = log[k]*PHI[k] / gcd and order expo / gcd
 * so that log = expo[k] */
void
acb_dirichlet_char_conrey(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t x)
{
    /* assume chi->x already set if x == NULL */
    if (x == NULL)
        x = chi->x;
    else
        acb_dirichlet_conrey_copy(chi->x, G, x);

    chi->q = G->q;
    chi->parity = acb_dirichlet_conrey_parity(G, x);
    chi->conductor = acb_dirichlet_conrey_conductor(G, x);

    acb_dirichlet_char_set_expo(chi, G);
    /* optional: divide by gcd to obtain true order */
    acb_dirichlet_char_normalize(chi, G);
}
