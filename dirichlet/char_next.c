/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "dirichlet.h"

int
dirichlet_fullchar_next(dirichlet_fullchar_t chi, const dirichlet_group_t G)
{
    int k;
    k = dirichlet_char_next(chi->x, G);
    dirichlet_fullchar_char(chi, G, NULL);
    /* return last index modified */
    return k;
}
