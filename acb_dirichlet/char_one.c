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
acb_dirichlet_char_one(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G)
{
    acb_dirichlet_conrey_one(chi->x, G);
    chi->q = G->q;
    chi->conductor = 1;
    chi->parity = 0;
    acb_dirichlet_char_set_expo(chi, G);
    acb_dirichlet_char_normalize(chi, G);
}
