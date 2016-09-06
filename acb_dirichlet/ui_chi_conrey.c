/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

ulong
acb_dirichlet_ui_chi_conrey(const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, const acb_dirichlet_conrey_t x)
{
        ulong v = 0, k;

        /* TODO: nmod_addmul? */
        for (k = 0; k < G->num; k++)
            v = nmod_add(v, nmod_mul(chi->expo[k], x->log[k], chi->order), chi->order);

        return v;
}
