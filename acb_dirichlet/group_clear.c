/*
    Copyright (C) 2015 Jonathan Bober
    Copyright (C) 2016 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

void
acb_dirichlet_group_clear(acb_dirichlet_group_t G)
{
    flint_free(G->primes);
    flint_free(G->exponents);
    flint_free(G->primepowers);
    flint_free(G->generators);
    flint_free(G->phi);
    flint_free(G->PHI);
}
