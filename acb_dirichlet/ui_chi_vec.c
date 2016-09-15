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
acb_dirichlet_ui_chi_vec(ulong *v, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong nv)
{
    if (2 * nv > G->q)
        acb_dirichlet_ui_chi_vec_loop(v, G, chi, nv);
    else
        acb_dirichlet_ui_chi_vec_primeloop(v, G, chi, nv);
}
