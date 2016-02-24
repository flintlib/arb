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

void
acb_dirichlet_char_init(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G) {
  chi->expo = flint_malloc(G->num * sizeof(ulong));
}

void
acb_dirichlet_char_clear(acb_dirichlet_char_t chi) {
    flint_free(chi->expo);
}

void
acb_dirichlet_char_log(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G, const acb_conrey_t x)
{
  ulong k, g;
  g = G->expo;
  /* modify exponents */
  for (k = 0; k < G->num; k++)
  {
    chi->expo[k] = (x->log[k] * G->PHI[k]) % G->expo;
    g = n_gcd(g, chi->expo[k]);
  }
  for (k = 0; k < G->num; k++)
    chi->expo[k] = chi->expo[k] / g;

  chi->q = G->q;
  chi->order = G->expo / g;
  chi->n = x->n;
}

void
acb_dirichlet_char(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G, ulong n)
{
  acb_conrey_t x;
  x->log = chi->expo;
  acb_conrey_log(x, G, n);
  acb_dirichlet_char_log(chi, G, x);
}

