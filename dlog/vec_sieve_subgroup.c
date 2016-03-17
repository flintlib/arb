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
#include <math.h>

#define vbs 0

void
dlog_vec_sieve_subgroup(ulong *v, ulong nv, ulong a, ulong va, ulong M, nmod_t mod, ulong na, nmod_t order)
{
    ulong smooth = 0, sievecount = 0, logcount = 0, missed = 0;
    ulong logcost, limcount;
    ulong k, p, p1, pmax, logm1;
    log_pair_t * w;
    dlog_precomp_t pre;
    n_primes_t iter;
    ulong X, aX, vaX;

    /* store size */
    w = flint_malloc( nv * sizeof(log_pair_t));
    for (k = 0; k < nv; k++)
    {
        w[k].m = 1;
        w[k].logm = 0; /* could be v[k]... */
    } 

    /* discrete log on first primes, then sieve */
    pmax = (nv < mod.n) ? nv : mod.n;
    p1 = 50; /* FIXME: tune this limit! */
    dlog_precomp_n_init(pre, a, mod.n, na, p1);
    /*flint_printf("## single log cost: %wu\n", pre->cost);*/
    logcost = pre->cost;

    if (logcost < 15)
    {
        /* p1 = pmax; */
        limcount = mod.n;
    }
    else
    {
        limcount = ceil(pow((double)mod.n,1./2.3) * 40 / logcost);
        logm1 = (mod.n % 2) ? 0 : dlog_precomp(pre, mod.n - 1);
    }

    /* find big power of gen */
    X = n_nextprime(na / 2 + 10, 0);
    X = (na % 257) ? 257 % na : 1031 % na ; /* FIXME! */
    aX = nmod_pow_ui(a, X, mod);
    vaX = nmod_mul(va, X % order.n, order);

    n_primes_init(iter);
    while ((p = n_primes_next(iter)) < pmax)
    { 
        double cost = log(mod.n)/log(p);
        ulong m, wp; 
        if (mod.n % p == 0) /* FIXME: those primes could be known... */
            continue; /* won't be attained another time */
        cost = log(mod.n)/log(p);
        cost = pow(cost,cost);
        sievecount++;
        /* if (p < p1 || (wp = logp_sieve(w, nv, p, mod.n, logm1, order, logcost)) == NOT_FOUND) */
        /*if (smooth < limcount || (wp = logp_sieve_factor(w, nv, p, mod.n, a, na, va, logm1, order, logcost)) == NOT_FOUND)*/
        if (logcost < cost || (wp = dlog_vec_pindex_factorgcd(w, nv, p, mod, aX, na, vaX, logm1, order, cost)) == NOT_FOUND)
        {
            if (logcost < cost)
                sievecount--;
            else
                missed++;
            logcount++;
            wp = nmod_mul(dlog_precomp(pre, nmod_pow_ui(p, M, mod)), va, order);
        }
        for (k = p, m = 1; k < nv; k += p, m++)
        {
            w[k].m = w[m].m * p;
            w[k].logm = nmod_add(w[m].logm,  wp, order);
            if (w[k].m == k)
                smooth++;
        }
    }
    /* write in v */
    for (k = 0; k < nv; k++)
        if (v[k] != NOT_FOUND)
            v[k] = nmod_add(v[k], w[k].logm, order);
#if vbs
    if (missed)
        flint_printf("[sieve: got %wu / %wu, n = %wu, cost %wu, logs %wu, sieve %wu missed %wu]\n",
                smooth, limcount, mod.n, logcost, logcount, sievecount, missed);
#endif
    n_primes_clear(iter);
    flint_free(w);
}
