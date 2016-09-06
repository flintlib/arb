/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

int
acb_dirichlet_char_eq(const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi1, const acb_dirichlet_char_t chi2)
{
    acb_dirichlet_conrey_t x, y;

    if (chi1->q != chi2->q)
        return 0;

    if (chi1->order.n != chi2->order.n)
        return 0;

    if (chi1->conductor != chi2->conductor)
        return 0;

    if (!acb_dirichlet_conrey_eq(G, chi1->x, chi2->x))
        return 0;

    x->n = y->n = 1;
    x->log = chi1->expo;
    y->log = chi2->expo;
    if (!acb_dirichlet_conrey_eq(G, x, y))
        return 0;

    return 1;
}
