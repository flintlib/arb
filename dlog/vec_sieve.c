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

#include "dlog.h"

static ulong
logp_sieve(log_pair_t * v, ulong nv, ulong p, ulong mod, ulong logm1, nmod_t order, int maxtry)
{
    int i, j, nr = 0, nm = 0, ng = 0;
    ulong l, pm, logm;
    ulong u[2], r[2], t;
#if rnd
    flint_rand_t state;
    flint_randinit(state);
#endif
#if vbs
    flint_printf("\nEnter logp_sieve p=%wu mod %wu...\n", p, mod);
#endif
    while (1) {
        /* find multiplier m */
        logm = 0;
        pm = l = p;
        do {
            nm++;
            /* random ? pb when p lies in a small subgroup */
            do {
                nr++;
#if rnd
                l = 1 + n_randint(state, p - 1);
#else
                l = (l > 1) ? l - 1 : p - 1;
#endif
            } while (v[l].m != l);
            pm *= l;
            logm += v[l].logm;
        } while (pm < mod);
        pm = pm % mod;
#if vbs
        flint_printf("[pm=%wu, v[pm]=%ld]", pm, v[pm]);
#endif
        /* half gcd u * pm + v * mod = r, ignore v */
        u[0] = 0; r[0] = mod;
        u[1] = 1; r[1] = pm;
        i = 1; j = 0; /* flip flap */
        do {
            ng++;
#if vbs
            flint_printf("[r=%d, v[r]=%d, u=%d, v[u]=%d]\n",
                    r[i],v[r[i]], u[i], v[u[i]]);
#endif
            if (r[i] < nv && v[r[i]].m == r[i] && u[i] < nv && v[u[i]].m == u[i])
            {
                ulong x;
                /* chi(-1)^j*chi(u)*chi(p)*chi(m)=chi(r) */
                x = nmod_sub(v[r[i]].logm, nmod_add(v[u[i]].logm, logm, order), order);
                if (j)
                    x = nmod_add(x, logm1, order);
#if rnd
                flint_randclear(state);
#endif
                return x;
            }

            j = i; i = 1 - i; /* switch */
            t = r[i] / r[j];
            r[i] = r[i] % r[j];
            u[i] = u[i] + t * u[j]; /* (-1)^j */

        } while (r[i] > 0 && u[i] < p);
        if (nm > maxtry)
            return NOT_FOUND;
    }
}

void
dlog_vec_sieve(ulong *v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order)
{
    int maxtry;
	ulong k, p, p1, pmax, logm1;
	log_pair_t * w;
	dlog_precomp_t pre;
	n_primes_t iter;

    /* store size */
    w = flint_malloc( nv * sizeof(log_pair_t));
	for (k = 0; k < nv; k++)
    {
	    w[k].m = 1;
	    w[k].logm = 1;
    }
	w[1].logm = 0; 

	/* discrete log on first primes, then sieve */
	pmax = (nv < mod.n) ? nv : mod.n;
	p1 = maxtry = 50; /* FIXME: tune this limit! */
	dlog_precomp_n_init(pre, a, mod.n, na, p1);

	logm1 = (mod.n % 2) ? 0 : dlog_precomp(pre, mod.n - 1);

	n_primes_init(iter);
	while ((p = n_primes_next(iter)) < pmax)
	{ 
        ulong m, wp; 
        if (mod.n % p == 0) /* FIXME: those primes could be known... */
            continue; /* won't be attained another time */
        if (p < p1 || (wp = logp_sieve(w, nv, p, mod.n, logm1, order, maxtry)) != NOT_FOUND)
            wp = nmod_mul(dlog_precomp(pre, p), va, order);
        for (k = p, m = 1; k < nv; k += p, m++)
        {
            w[k].m *= p;
            w[k].logm = nmod_add(w[k].logm,  wp, order);
        }
    }
    /* write in v */
	for (k = 0; k < nv; k++)
        if (v[k] != NOT_FOUND)
            v[k] = nmod_add(v[k], w[k].logm, order);
    n_primes_clear(iter);
    flint_free(w);
}
