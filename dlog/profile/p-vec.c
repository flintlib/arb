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
#include "profiler.h"
typedef void (*vec_f) (ulong *v, ulong nv, ulong a, ulong va, const nmod_t mod, ulong na, const nmod_t order);
   
void
f_empty(ulong *v, ulong nv, ulong a, ulong va, const nmod_t mod, ulong na, const nmod_t order)
{
    return;
}

int main()
{
    int i, ni = 3;
    int bits[5] = { 10, 20, 30, 40, 50 };

    int j, nj = 5;
    ulong * v;
    ulong nv[5] = { 50, 200, 1000, 2000, 10000 };

    int k, np = 1000;
    nmod_t * p;
    ulong * a;

    int l, nf = 3;
    vec_f func[4] = { f_empty, dlog_vec_loop, dlog_vec_sieve, dlog_vec_crt };
    char * n[4] = { "empty", "loop", "sieve", "crt" };

    flint_rand_t state;
    nmod_t order;

    nmod_init(&order, 100);
    p = flint_malloc(np * sizeof(nmod_t));
    a = flint_malloc(np * sizeof(ulong));

    for (i = 0; i < ni; i++)
    {
        for (k = 0; k < np; k++)
        {
            nmod_init(&p[k], n_randprime(state, bits[i], 0));
            a[k] = n_primitive_root_prime(p[k].n);
        }

        for (j = 0; j < nj; j++)
        {

            v = flint_malloc(nv[j] * sizeof(ulong));

            flint_printf("log(1..%wu) mod %d primes of size %d bits....\n", nv[j], np, bits[i]);
            fflush(stdout);

            for (l = 0; l < nf; l++)
            { 
                if (l == 1 && (i >= 2 || j > 0))
                    continue;
                
                flint_printf("%-20s...   ",n[l]);
                fflush(stdout);
                TIMEIT_ONCE_START
                for (k = 0; k < np; k++)
                {
                    int kk;
                    for (kk=0; kk < nv[j]; kk++)
                        v[kk] = 0;
                    (func[l])(v, nv[j], a[k], 1, p[k], p[k].n - 1, order);
                }
                TIMEIT_ONCE_STOP
            }
            flint_free(v);
        }
    }

    flint_free(p);
    flint_free(a);
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
