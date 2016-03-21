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

typedef void (*vec_f) (ulong *v, ulong nv, ulong a, ulong va, const nmod_t mod, ulong na, const nmod_t order);
   
void
dlog_vec_trivial(ulong *v, ulong nv, ulong a, ulong va, const nmod_t mod, ulong na, const nmod_t order)
{
    ulong k;
    dlog_precomp_t pre;
    dlog_precomp_n_init(pre, a, mod.n, na, 50);
    for (k = 1; k < nv; k++)
        if (n_gcd(k, mod.n) > 1)
            v[k] = DLOG_NOT_FOUND;
        else
            v[k] = dlog_precomp(pre, k % mod.n);
    dlog_precomp_clear(pre);
}

static ulong
dlog_vec_diff(ulong * v, ulong * ref, ulong nv)
{
    ulong k;
    for (k = 1; k < nv; k++)
        if (ref[k] != v[k])
            return k;
    return 0;
}

int main()
{
    slong bits, nv, iter;
    flint_rand_t state;
    int f, nf = 4;
    vec_f func[4] = { dlog_vec_trivial, dlog_vec_loop, dlog_vec_eratos,
        dlog_vec_sieve };
    char * n[4] = { "trivial", "loop", "eratos", "sieve" };


    flint_printf("dlog(1..n) mod p....");
    fflush(stdout);
    flint_randinit(state);

    for (bits = 10; bits <= 35; bits += 5)
    {
        for (nv = 10; nv <= 10000; nv *= 10)
        {

            ulong *v, *ref;
            int iref;

            iref = (bits == 10 && nv <= 1000) ? 0 : 2;

            ref = flint_malloc(nv * sizeof(ulong));
            v = flint_malloc(nv * sizeof(ulong));

            for (iter = 0; iter < 10; iter++)
            {

                int k;
                ulong p, a, va, na;
                nmod_t mod, order;

                p = n_randprime(state, bits, 0);
                a = n_primitive_root_prime(p);

                nmod_init(&mod, p);
                va = 1; na = p - 1;
                nmod_init(&order, na);

                dlog_vec_fill(ref, nv, 0);
                (func[iref])(ref, nv, a, va, mod, na, order);

                /* compare */
                for (f = iref + 1; f < nf; f++)
                {
                    dlog_vec_fill(v, nv, 0);
                    (func[f])(v, nv, a, va, mod, na, order);

                    if ((k = dlog_vec_diff(v, ref, nv)))
                    {
                        flint_printf("FAIL: log(%wu,%wu) mod %wu: %s->%w != %s->%w\n",
                                k, a, p, n[iref], ref[k], n[f], v[k]);
                        abort();
                    }
                }
            }
            flint_free(ref);
            flint_free(v);
        }
    }
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
