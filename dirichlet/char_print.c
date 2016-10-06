/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "dirichlet.h"

void
dirichlet_char_print(const dirichlet_group_t G, const dirichlet_char_t chi)
{
    dirichlet_conrey_t x;
    flint_printf("chi_%wu(%wu,.) of order %wu, parity %wd, index ", G->q, chi->x->n, chi->order, chi->parity);
    dirichlet_conrey_print(G, chi->x);
    flint_printf(" and exponents ");
    x->log = chi->expo;
    dirichlet_conrey_print(G, x);
}
