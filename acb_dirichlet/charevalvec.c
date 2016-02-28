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

static void
set_non_invertible_values(long *v, const acb_dirichlet_group_t G, ulong nv)
{
  ulong k, l;
  if (G->q_even > 1)
  {
    for (k = 2; k < nv; k += 2)
      v[k] = -1;
  }
  for (l = 0; l < G->num; l++)
  {
    ulong p = G->primes[k];
    for (k = p; k < nv; k += p)
      v[k] = -1;
  }
}

/* loop over whole group */
void
n_dirichlet_char_vec_large(long *v, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, ulong nv)
{
  ulong t, k, j;
  acb_conrey_t x;
  acb_conrey_init(x, G);
  acb_conrey_one(x, G);
  t = v[1] = 0;
  while ( (j = acb_conrey_next(x, G)) < G->num )
  {
    /* exponents were modified up to j */
    for (k = 0; k < j; k++)
      t = (t + chi->expo[k] * x->log[k]) % chi->order;
    if (x->n < nv)
      v[x->n] = t;
  }
  /* fix result outside primes */
  set_non_invertible_values(v, G, nv);
  /* copy outside modulus */
  for (k = G->q + 1; k < nv ; k++ )
    v[k] = v[k-G->q];
  acb_conrey_clear(x);
}

/* loop over primary components */
void
n_dirichlet_char_vec_med(long *v, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, ulong nv)
{
  /* prime by prime */
  ulong k, l;

  for(k = 1; k < nv; ++k)
      v[k] = 0;

  for(l = 1; l < G->num; ++l)
  {
    long p, pe, g, x, vp, xp;
    long j, vj;
    p = G->primes[l];
    pe = G->primepowers[l];
    g = G->generators[l] % pe;
    vj = vp = chi->expo[l];
    if( g < 0 ) g += pe;
    /* for each x = g^j mod p^e,
     * set a[x] += j*vp
     * and use periodicity */
    for(j = 1, x = g; x > 1; j++)
    {
      for(xp = x; xp < nv; xp+=pe)
          v[xp] = (v[xp] + vj) % chi->order;
      x = (x*g) % pe;
      vj = (vj + vp) % chi->order;
    }
  }
  set_non_invertible_values(v, G, nv);
}


/* eratosthene sieve */
void
n_dirichlet_char_vec_small(long *v, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, ulong nv)
{
	ulong k, p, pmax;
	n_primes_t iter;

	n_primes_init(iter);

	pmax = (nv < G->q) ? nv : G->q;
	v[1] = 0; 
	while ((p = n_primes_next(iter)) < pmax)
	{ 
		if (G->q % p == 0) 
		{
			for (k = p; k < nv; k += p)
				v[k] = -1;
		}
		else
		{ 
			long chip; 
			chip = n_dirichlet_char_eval(G, chi, p);
			for (k = p; k < nv; k += p)
				if (v[k] != -1)
					 v[k] = (v[k] + chip) % chi->order;
		}
	}
	n_primes_clear(iter);
}
