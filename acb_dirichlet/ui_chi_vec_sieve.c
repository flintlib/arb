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

    Copyright (C) 2016 Pascal Molin

******************************************************************************/

#include "acb_dirichlet.h"

/* sieve on primes */
/* TODO: see if really more efficient than primeloop using dlog
 * this one is cheaper propagating the values, but not sure this
 * is noticeable ...
 */
void
acb_dirichlet_ui_chi_vec_sieve(ulong *v, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong nv)
{
	slong k, p, pmax;
	n_primes_t iter;

	n_primes_init(iter);

	pmax = (nv < G->q) ? nv : G->q;
	v[1] = 0; 

	while ((p = n_primes_next(iter)) < pmax)
	{ 
		if (G->q % p == 0) 
		{
			for (k = p; k < nv; k += p)
				v[k] = ACB_DIRICHLET_CHI_NULL;
		}
		else
		{ 
			long chip; 
			chip = acb_dirichlet_ui_chi(G, chi, p);

			for (k = p; k < nv; k += p)
				if (v[k] != -1)
					 v[k] = nmod_add(v[k], chip, chi->order);
		}
	}

	n_primes_clear(iter);
}
