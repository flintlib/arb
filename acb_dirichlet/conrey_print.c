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
acb_dirichlet_conrey_print(const acb_dirichlet_group_t G, const acb_dirichlet_conrey_t x)
{
    slong k;
    if (G->num)
        flint_printf("[%wu", x->log[0]);
    else
        flint_printf("[");
    for (k = 1; k < G->num; k++)
        flint_printf(", %wu", x->log[k]);
    flint_printf("]");
}
