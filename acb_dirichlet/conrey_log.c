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

#include "dlog.h"
#include "acb_dirichlet.h"

/* TODO: use dlog module instead of n_discrete_log_bsgs */
/* assume m is invertible */
void
acb_dirichlet_conrey_log(acb_dirichlet_conrey_t x, const acb_dirichlet_group_t G, ulong m)
{
    ulong k, pk, gk;
    /* even part */
    if (G->neven >= 1)
    {
      x->log[0] = (m % 4 == 3);
      if (G->neven == 2)
      {
        ulong m2 = (x->log[0]) ? -m % G->q_even : m % G->q_even;
        if (G->P[1].dlog == NULL)
            x->log[1] = dlog_mod2e_1mod4(m2, G->P[1].e,
                    nmod_inv(5, G->P[1].pe), G->P[1].pe);
        else
            x->log[1] = dlog_precomp(G->P[1].dlog, m2);
      }
    }
    /* odd part */
    for (k = G->neven; k < G->num; k++)
    {
        if (G->P[k].dlog == NULL)
        {
            pk = G->P[k].pe.n;
            gk = G->P[k].g;
            x->log[k] = n_discrete_log_bsgs(m % pk, gk, pk);
        }
        else
        {
            x->log[k] = dlog_precomp(G->P[k].dlog, m % G->P[k].pe.n);
        }
    }
    /* keep value m */
    x->n = m;
}
