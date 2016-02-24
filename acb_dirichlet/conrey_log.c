/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2015 Jonathan Bober
    Copyright (C) 2016 Fredrik Johansson

******************************************************************************/

#include "acb_dirichlet.h"

/* todo: modular arithmetic
   discrete log can be computed along exponents or using p-adic log
*/

void
acb_conrey_init(acb_conrey_t x, const acb_dirichlet_group_t G) {
    x->log = flint_malloc(G->num * sizeof(ulong));
}

void
acb_conrey_one(acb_conrey_t x, const acb_dirichlet_group_t G) {
    ulong k;
    for (k = 0; k < G->num ; k++)
      x->log[k] = 0;
    x->n = 1;
}

void
acb_conrey_clear(acb_conrey_t x) {
    flint_free(x->log);
}

/* TODO: use precomputations in G if present */
void
acb_conrey_log(acb_conrey_t x, const acb_dirichlet_group_t G, ulong m)
{
    ulong k, pk, gk;
    /* even part */
    if (G->neven >= 1)
      x->log[0] = (m % 4 == 3);
    if (G->neven == 2)
    {
        ulong q_even = G->q_even;
        ulong g2 = 5;
        ulong m2 = (m % 4 == 3) ? n_negmod(m, q_even) : m % q_even;
        x->log[1] = n_discrete_log_bsgs(m2, g2, q_even);
    }
    /* odd part */
    for (k = G->neven; k < G->num; k++)
    {
        pk = n_pow(G->primes[k], G->exponents[k]);
        gk = G->generators[k] % pk;
        x->log[k] = n_discrete_log_bsgs(m % pk, gk, pk);
    }
    /* keep value m */
    x->n = m;
}

int
acb_conrey_next(acb_conrey_t x, const acb_dirichlet_group_t G)
{
  /* update index */
  ulong k;
  for (k=0; k < G->num ; k++)
  {
    x->n = x->n * G->generators[k];
    if(x->log[k]++ < G->phi[k])
      break;
    x->log[k] = 0;  
  }
  /* return last index modified */
  return k;
}
