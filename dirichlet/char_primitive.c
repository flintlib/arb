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
dirichlet_char_primitive(dirichlet_char_t chi0, const dirichlet_group_t G0, const dirichlet_group_t G, const dirichlet_char_t chi)
{
    chi0->q = chi->conductor;
    chi0->parity = chi->parity;
    chi0->conductor = chi->conductor;
    dirichlet_conrey_primitive(chi0->x, G, chi->x, chi->conductor);
    dirichlet_char_set_expo(chi0, G0);
    /* optional: divide by gcd to obtain true order */
    dirichlet_char_normalize(chi0, G0);
}
