/*
    Copyright (C) 2015 Jonathan Bober
    Copyright (C) 2016 Fredrik Johansson
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "dirichlet.h"

/* todo: modular arithmetic */

ulong
dirichlet_ui_pairing_conrey(const dirichlet_group_t G, const dirichlet_conrey_t a, const dirichlet_conrey_t b)
{
    ulong x, k;
    x = 0;

    for (k = 0; k < G->num; k++)
        x = n_addmod(x, G->PHI[k] * n_mulmod2(a->log[k], b->log[k], G->P[k].phi), G->expo);

    return x;
}
