/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

void
acb_dirichlet_conrey_first_primitive(acb_dirichlet_conrey_t x, const acb_dirichlet_group_t G)
{
    ulong k;
    if (G->q % 4 == 2)
    {
        flint_printf("Exception (acb_dirichlet_conrey_first_primitive). No primitive element mod %wu.\n",G->q);
        abort();
    }
    x->n = 1;
    for (k = 0; k < G->num ; k++)
    {
        if (k == 0 && G->neven == 2)
            x->log[k] = 0;
        else
        {
            x->n = nmod_mul(x->n, G->generators[k], G->mod);
            x->log[k] = 1;
        }
    }
}
