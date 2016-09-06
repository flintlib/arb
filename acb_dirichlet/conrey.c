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
acb_dirichlet_conrey_init(acb_dirichlet_conrey_t x, const acb_dirichlet_group_t G) {
    x->log = flint_malloc(G->num * sizeof(ulong));
}

void
acb_dirichlet_conrey_clear(acb_dirichlet_conrey_t x) {
    flint_free(x->log);
}
