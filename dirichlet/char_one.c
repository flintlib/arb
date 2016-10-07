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
dirichlet_fullchar_one(dirichlet_fullchar_t chi, const dirichlet_group_t G)
{
    slong k;
    dirichlet_char_one(chi->x, G);
    chi->q = G->q;
    chi->conductor = 1;
    chi->parity = 0;
    for (k = 0; k < G->num; k++)
        chi->expo[k] = 0;
    nmod_init(&(chi->order), 1);
}
