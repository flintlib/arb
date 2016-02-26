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

/* FIXME: multiplications mod G->q should be n_mulmod to avoid overflow */

void
acb_dirichlet_char_init(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G) {
  chi->expo = flint_malloc(G->num * sizeof(ulong));
}

void
acb_dirichlet_char_clear(acb_dirichlet_char_t chi) {
    flint_free(chi->expo);
}

static void
acb_dirichlet_char_normalize(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G)
{
  ulong k, g;
  g = G->expo;
  for (k = 0; k < G->num; k++)
    g = n_gcd(g, chi->expo[k]);
  for (k = 0; k < G->num; k++)
    chi->expo[k] = chi->expo[k] / g;
  chi->order = G->expo / g;
}

static void
acb_dirichlet_char_denormalize(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G)
{
  ulong k, g;
  g = G->expo / chi->order;
  for (k = 0; k < G->num; k++)
    chi->expo[k] *= g;
}

/* char n has exponents  = log[k]*PHI[k] / gcd and order expo / gcd 
 * so that log = expo[k] */
static void
acb_dirichlet_char_conrey(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G, const acb_conrey_t x)
{
  ulong k;
  chi->q = G->q;
  chi->n = x->n;

  for (k = 0; k < G->num; k++)
    chi->expo[k] = (x->log[k] * G->PHI[k]) % G->expo;

  /* optional: divide by gcd to obtain true order */
  acb_dirichlet_char_normalize(chi, G);
}

void
acb_dirichlet_char(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G, ulong n)
{
  acb_conrey_t x;
  x->log = chi->expo;
  acb_conrey_log(x, G, n);
  acb_dirichlet_char_conrey(chi, G, x);
}

ulong
acb_dirichlet_char_next(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G)
{
  ulong k;

  acb_dirichlet_char_denormalize(chi, G);

  /* update index */
  for (k=0; k < G->num ; k++)
  {
    /* chi->n = n_mulmod(chi->n, G->generators[k], G->q); */
    chi->n = chi->n * G->generators[k] % G->q;
    chi->expo[k] += G->PHI[k];
    if (chi->expo[k] < G->expo)
      break;
    chi->expo[k] = 0;  
  }

  acb_dirichlet_char_normalize(chi, G);

  /* return last index modified */
  return k;
}

ulong
acb_dirichlet_char_next_primitive(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G)
{
  ulong k;

  acb_dirichlet_char_denormalize(chi, G);

  /* update index */
  k = 0;
  if (G->neven == 2)
  {
      /* chi->n = n_mulmod(chi->n, G->generators[0], G->q); */
      chi->n = chi->n * G->generators[0] % G->q;
      if (++chi->expo[0] < G->expo)
          return 0;
      chi->expo[0] = 0;
      k = 1;
  }
  for (; k < G->num ; k++)
  {
      /* chi->n = n_mulmod(chi->n, G->generators[k], G->q); */
      chi->n = chi->n * G->generators[k] % G->q;
      chi->expo[k] += G->PHI[k];
      if (chi->expo[k] % G->primes[k] == 0)
      {
          /* chi->n = n_mulmod(chi->n, G->generators[k], G->q); */
          chi->n = chi->n * G->generators[k] % G->q;
          chi->expo[k] += G->PHI[k];
      }
      if (chi->expo[k] < G->expo)
          break;
      chi->expo[k] = G->PHI[k];  
  }

  acb_dirichlet_char_normalize(chi, G);

  /* return last index modified */
  return k;
}

void
acb_dirichlet_char_one(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G)
{
  ulong k;
  chi->q = G->q;
  chi->n = 1;
  for (k = 0; k < G->num; k++)
      chi->expo[k] = 0;
  chi->order = 1;
}

void
acb_dirichlet_char_first_primitive(acb_dirichlet_char_t chi, const acb_dirichlet_group_t G)
{
  acb_conrey_t x;
  chi->q = G->q;
  x->log = chi->expo;
  acb_conrey_first_primitive(x, G);
  chi->n = x->n;
  acb_dirichlet_char_normalize(chi, G);
}
