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

#define vbs 1

void
dlog_vec_sieve(ulong *v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order)
{
    ulong smooth = 0, sievecount = 0, logcount = 0, missed = 0;
    ulong logcost, limcount;
    ulong k, p, p1, pmax, logm1;
    dlog_precomp_t pre;
    n_primes_t iter;
    ulong X, aX, vaX;

    dlog_vec_fill(v, nv, DLOG_NOT_FOUND);
    v[1] = 0;

    logm1 = (na % 2) ? 0 : nmod_mul(na / 2, va, order);

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
        ulong m, vp; 
        if (mod.n % p == 0)
            continue; /* won't be attained another time */
        cost = log(mod.n)/log(p);
        cost = pow(cost,cost);
        sievecount++;
        /* if (p < p1 || (wp = logp_sieve(w, nv, p, mod.n, logm1, order, logcost)) == NOT_FOUND) */
        /*if (smooth < limcount || (wp = logp_sieve_factor(w, nv, p, mod.n, a, na, va, logm1, order, logcost)) == NOT_FOUND)*/
        if (logcost < cost || (vp = dlog_vec_pindex_factorgcd(v, nv, p, mod, aX, na, vaX, logm1, order, cost)) == DLOG_NOT_FOUND)
        {
            if (logcost < cost)
                sievecount--;
            else
                missed++;
            logcount++;
            vp = nmod_mul(dlog_precomp(pre, p), va, order);
        }
        for (k = p, m = 1; k < nv; k += p, m++)
        {
            if (v[m] == DLOG_NOT_FOUND)
                continue;
            smooth++;
            v[k] = nmod_add(v[m],  vp, order);
        }
    }
#if vbs
    if (missed)
        flint_printf("[sieve: got %wu / %wu, n = %wu, cost %wu, logs %wu, sieve %wu missed %wu]\n",
                smooth, limcount, mod.n, logcost, logcount, sievecount, missed);
#endif
    n_primes_clear(iter);
    for (k = mod.n + 1; k < nv; k++)
        v[k] = v[k - mod.n];
}
