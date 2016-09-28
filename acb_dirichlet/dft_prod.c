/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

acb_dirichlet_dft_step_ptr
_acb_dirichlet_dft_steps_prod(slong * cyc, slong num, slong prec)
{
    slong i, len, dv;
    acb_dirichlet_dft_step_ptr s;
    s = flint_malloc(num * sizeof(acb_dirichlet_dft_step_struct));

    len = 1; dv = 1;
    for (i = 0; i < num; i++)
        len *= cyc[i];

    for (i = 0; i < num; i++)
    {
        slong m = cyc[i];
        len /= m;
        s[i].m = m;
        s[i].M = len;
        s[i].dv = len;
        s[i].dz = 0;
        s[i].z = NULL;
        _acb_dirichlet_dft_precomp_init(s[i].pre, m, NULL, 0, len, prec);
        dv *= m;
    }

    return s;
}

void
acb_dirichlet_dft_prod_clear(acb_dirichlet_dft_prod_t t)
{
    slong i;
    for (i = 0; i < t->num; i++)
        acb_dirichlet_dft_precomp_clear(t->cyc[i].pre);
    flint_free(t->cyc);
}


void
acb_dirichlet_dft_prod_precomp(acb_ptr w, acb_srcptr v, const acb_dirichlet_dft_prod_t prod, slong prec)
{
    if (prod->num == 0)
        acb_set(w + 0, v + 0);
    else
        acb_dirichlet_dft_step(w, v, prod->cyc, prod->num, prec);
}

void
acb_dirichlet_dft_prod(acb_ptr w, acb_srcptr v, slong * cyc, slong num, slong prec)
{
    acb_dirichlet_dft_prod_t t;
    acb_dirichlet_dft_prod_init(t, cyc, num, prec);
    acb_dirichlet_dft_prod_precomp(w, v, t, prec);
    acb_dirichlet_dft_prod_clear(t);
}
