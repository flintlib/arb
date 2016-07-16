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
#include "profiler.h"
typedef void (*dir_f) (ulong *v, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong nv);

static void
dir_empty(ulong *v, const acb_dirichlet_group_t G, const acb_dirichlet_char_t chi, slong nv)
{
    return;
}

static void
vecloop(dir_f dir, ulong minq, ulong maxq, ulong * rand, ulong nr, ulong * v, ulong nv)
{
    ulong q;
    TIMEIT_ONCE_START

    for (q = minq; q <= maxq; q++)
    {
        ulong r;
        acb_dirichlet_group_t G;
        acb_dirichlet_char_t chi;

        acb_dirichlet_group_init(G, q);
        acb_dirichlet_char_init(chi, G);

        for (r = 0; r < nr; r++)
        {
            acb_dirichlet_char(chi, G, rand[r] % q);
            dir(v, G, chi, nv);
        }

        acb_dirichlet_char_clear(chi);
        acb_dirichlet_group_clear(G);
    }

    TIMEIT_ONCE_STOP
    flint_printf("\n");
}

static void
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

int main()
{
    slong iter, k, nv, nref, r, nr;
    ulong minq, maxq;
    ulong * rand;
    int i, ni = 5;
    ulong q[5] =  {   2, 1000, 3000, 10000, 100000 };
    ulong qq[5] = { 500, 2000, 5000, 12000, 100500 };
    ulong * v;
    flint_rand_t state;

    nr = 20;

    flint_randinit(state);
    rand = flint_malloc(nr * sizeof(ulong));
    v = flint_malloc(nv * sizeof(ulong));

    for (r = 0; r < nr; r++)
        rand[r] = n_randprime(state, 42, 0);

    for (i = 0; i < ni; i++)
    {

        ulong minq = q[i], maxq = qq[i];
        nv = 2000;

        flint_printf("%wu * chi(rand, 1..%wu) for all %wu <= q <= %wu....\n", nr, nv, minq, maxq);
        fflush(stdout);

        flint_printf("character only.......... ");
        fflush(stdout);
        vecloop(dir_empty, minq, maxq, rand, nr, v, nv);

        flint_printf("big loop................ ");
        fflush(stdout);
        vecloop(acb_dirichlet_ui_chi_vec_loop, minq, maxq, rand, nr, v, nv);

        flint_printf("med loop................ ");
        fflush(stdout);
        vecloop(acb_dirichlet_ui_chi_vec_primeloop, minq, maxq, rand, nr, v, nv);

        flint_printf("sieve................... ");
        fflush(stdout);
        vecloop(acb_dirichlet_ui_chi_vec_sieve, minq, maxq, rand, nr, v, nv);

        flint_printf("generic................. ");
        fflush(stdout);
        vecloop(acb_dirichlet_ui_chi_vec, minq, maxq, rand, nr, v, nv);
    }

    flint_free(v);
    flint_free(rand);
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
