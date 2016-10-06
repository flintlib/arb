/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "dirichlet.h"

ulong
dirichlet_ui_pairing(const dirichlet_group_t G, ulong m, ulong n)
{
    ulong x;
    dirichlet_conrey_t a, b;

    if (n_gcd(G->q, m) > 1 || n_gcd(G->q, n) > 1)
        return DIRICHLET_CHI_NULL;

    dirichlet_conrey_init(a, G);
    dirichlet_conrey_init(b, G);
    dirichlet_conrey_log(a, G, m);
    dirichlet_conrey_log(b, G, n);

    x = dirichlet_ui_pairing_conrey(G, a, b);

    dirichlet_conrey_clear(a);
    dirichlet_conrey_clear(b);

    return x;
}
