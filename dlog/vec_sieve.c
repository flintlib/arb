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
#define FACTOR_RATIO 4

static ulong
logp_sieve_gcd(log_pair_t * v, ulong nv, ulong p, nmod_t mod, ulong a, ulong na, ulong loga, ulong logm1, nmod_t order, int maxtry)
{
    int i, j, nm = 0, ng = 0;
    ulong pm, logm;
    ulong u[2], r[2], t;
#if vbs > 1
    flint_printf("\nEnter logp_sieve p=%wu mod %wu...\n", p, mod);
#endif
    pm = p;
    logm = 0;
    while (nm++ < maxtry)
    {
        pm = nmod_mul(pm, a, mod);
        logm = nmod_add(logm, loga, order);
        if (2 * pm > mod.n)
        {
            pm = nmod_neg(pm, mod);
            logm = nmod_add(logm, logm1, order);
        }
        /* half gcd u * pm + v * mod = r, ignore v */
        u[0] = 0; r[0] = mod.n;
        u[1] = 1; r[1] = pm;
        i = 1; j = 0; /* flip flap */
        do {
            ng++;
#if vbs > 1
            if (r[i] < nv && u[i] < nv)
                flint_printf("[r=%wu, v[r]=%wu, u=%wu, v[u]=%wu -- nv=%wu]\n",
                        r[i], v[r[i]].m, u[i], v[u[i]].m, nv);
#endif
            if (r[i] < nv && v[r[i]].m == r[i] && u[i] < nv && v[u[i]].m == u[i])
            {
                ulong x;
                /* chi(-1)^j*chi(u)*chi(p)*chi(m)=chi(r) */
                x = nmod_sub(v[r[i]].logm, nmod_add(v[u[i]].logm, logm, order), order);
                if (j)
                    x = nmod_add(x, logm1, order);
                return x;
            }

            j = i; i = 1 - i; /* switch */
            t = r[i] / r[j];
            r[i] = r[i] % r[j];
            u[i] = u[i] + t * u[j]; /* (-1)^j */

        } while (r[i] > 0 && u[i] < nv);
    }
    return NOT_FOUND;
}

static ulong
logp_sieve_factor(log_pair_t * v, ulong nv, ulong p, nmod_t mod, ulong a, ulong na, ulong loga, ulong logm1, nmod_t order, int maxtry)
{
    int nm = 0;
    ulong pm, logm;

    const ulong * prime;
    prime = n_primes_arr_readonly(p);

    pm = p;
    logm = 0;
    while (pm < mod.n)
    {
        /*nm++;*/ /* init ignored */
        pm *= a;
        logm = nmod_add(logm, loga, order);
    }
    pm = pm % mod.n;
    do {
        int i, j, ind[15], exp[15];
        /* find multiplier m */
        if (2 * pm > mod.n)
        {
            pm = nmod_neg(pm, mod);
            logm = nmod_add(logm, logm1, order);
        }
        
        for (i = 0, j = 0; j < p && pm >= nv && 4 * prime[j] < p; j++)
        {
            int e = n_remove(&pm, prime[j]);
            if (e)
            {
                ind[i] = j;
                exp[i] = e;
                i++;
            }
        }
        if (pm < nv && v[pm].m == pm)
        {
            /* goal! */
            ulong x = v[pm].logm;
            /* chi(m)*chi(p)=chi(pm)*prod */
            for(j = 0; j < i; j++)
                x = nmod_add(x, nmod_mul(exp[j], v[prime[ind[j]]].logm, order), order);
            x = nmod_sub(logm, x, order);
            /*flint_printf("managed %d mults / %d for p=%wu, pm=%wu\n",nm,maxtry,p,pm);*/
            return x;
        }
        nm++;
        pm = nmod_mul(pm, a, mod);
        logm = nmod_add(logm, loga, order);
    } while (nm < maxtry);
    return NOT_FOUND;
}

static int
factor_until(ulong * n, ulong nlim, const ulong * p, ulong pmax, ulong * fp, int * fe)
{
    int i, j;
    for (i = 0, j = 0; *n >= nlim && p[j] < pmax; j++)
        {
            int e = n_remove(n, p[j]);
            if (e)
            {
                fp[i] = p[j];
                fe[i] = e;
                i++;
            }
        }
    return i;
}

static ulong
logp_sieve_factorgcd(log_pair_t * v, ulong nv, ulong p, nmod_t mod, ulong a, ulong na, ulong loga, ulong logm1, nmod_t order, int maxtry)
{
    int nm = 0, ng = 0;
    ulong pm, logm, pmax;
    ulong u[2], r[2], t;
    ulong up[15], rp[15];
    int ue[15], re[15];
    const ulong * prime;
    prime = n_primes_arr_readonly(p);
    pmax = p / FACTOR_RATIO;
    pm = p;
    logm = 0;
    while (nm++ < maxtry)
    {
        int i, j, iu, ir;
        ulong logr;
        pm = nmod_mul(pm, a, mod);
        logm = nmod_add(logm, loga, order);
        /*
        if (2 * pm > mod.n)
        {
            pm = nmod_neg(pm, mod);
            logm = nmod_add(logm, logm1, order);
        }
        */
        /* half gcd u * pm + v * mod = r, ignore v */
        u[0] = 0; r[0] = mod.n;
        u[1] = 1; r[1] = pm;
        i = 1; j = 0; /* flip flap */
        while (r[i] > u[i])
        {
            ng++;
            if (r[i] < nv && v[r[i]].m == r[i] && u[i] < nv && v[u[i]].m == u[i])
            {
                /* early smooth detection: occurs for primes < 30 bits */
                ulong x;
                /* chi(-1)^j*chi(u)*chi(p)*chi(m)=chi(r) */
                x = nmod_sub(v[r[i]].logm, nmod_add(v[u[i]].logm, logm, order), order);
                if (j)
                    x = nmod_add(x, logm1, order);
                return x;
            }
            j = i; i = 1 - i; /* switch */
            t = r[i] / r[j];
            r[i] = r[i] % r[j];
            u[i] = u[i] + t * u[j]; /* (-1)^j */

        };
        /* try to factor both r[i] and u[i] */
        iu = factor_until(&u[i], nv, prime, pmax, up, ue);
        if (u[i] >= nv || v[u[i]].m < u[i])
            continue;
        ir = factor_until(&r[i], nv, prime, pmax, rp, re);
        if (r[i] >= nv || v[r[i]].m < r[i])
            continue;
        /* log(u)+log(p)+log(m)=log(r) */
        logr = 0;
        for (i=0; i < ir; i++)
            logr = nmod_add(logr, (re[i] * v[rp[i]].logm) % order.n, order);
        for (i=0; i < iu; i++)
            logm = nmod_add(logm, (ue[i] * v[up[i]].logm) % order.n, order);
        return nmod_sub(logr, logm, order);
    }
    return NOT_FOUND;
}

void
dlog_vec_sieve_ph(ulong *v, ulong nv, ulong a, ulong va, ulong M, nmod_t mod, ulong na, nmod_t order)
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
        if (logcost < cost || (wp = logp_sieve_factorgcd(w, nv, p, mod, aX, na, vaX, logm1, order, cost)) == NOT_FOUND)
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

void
dlog_vec_sieve(ulong *v, ulong nv, ulong a, ulong va, nmod_t mod, ulong na, nmod_t order)
{
    dlog_vec_sieve_ph(v, nv, a, va, 1, mod, na, order);
}
