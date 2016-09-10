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
acb_dirichlet_ui_chi(const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, ulong n)
{
    if (n_gcd(G->q, n) > 1)
    {
        return ACB_DIRICHLET_CHI_NULL;
    }
    else
    {
        ulong v;
        acb_dirichlet_conrey_t x;
        acb_dirichlet_conrey_init(x, G);
        acb_dirichlet_conrey_log(x, G, n);

        v = acb_dirichlet_ui_chi_conrey(G, chi, x);

        acb_dirichlet_conrey_clear(x);
        return v;
    }
}
