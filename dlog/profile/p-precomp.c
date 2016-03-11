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
#include <math.h>

typedef void (*multilog_f) (ulong p, ulong a, ulong num);

void
flog_table(ulong p, ulong a, ulong num)
{
    int k;
    dlog_table_t t;
    dlog_table_init(t, a, p);
    for (k = 1; k < num; k++)
    {
        if (k == p) continue;
        dlog_table(t, k % p);
    }
    dlog_table_clear(t);
}
void
flog_bsgs(ulong p, ulong a, ulong num)
{
    int k;
    dlog_bsgs_t t;
    dlog_bsgs_init(t, a, p, p-1, ceil( sqrt((double)num * p)));
    for (k = 1; k < num; k++)
    {
        if (k == p) continue;
        dlog_bsgs(t, k % p);
    }
    dlog_bsgs_clear(t);
}
void
flog_rho(ulong p, ulong a, ulong num)
{
    int k;
    dlog_rho_t t;
    dlog_rho_init(t, a, p, p-1);
    for (k = 1; k < num; k++)
    {
        if (k == p) continue;
        dlog_rho(t, k % p);
    }
    dlog_rho_clear(t);
}
void
flog_crt(ulong p, ulong a, ulong num)
{
    int k;
    dlog_crt_t t;
    dlog_crt_init(t, a, p, p-1, num);
    for (k = 1; k < num; k++)
    {
        if (k == p) continue;
        dlog_crt(t, k % p);
    }
    dlog_crt_clear(t);
}
void
flog_gen(ulong p, ulong a, ulong num)
{
    int k;
    dlog_precomp_t t;
    dlog_precomp_n_init(t, a, p, p-1, num);
    for (k = 1; k < num; k++)
    {
        if (k == p) continue;
        dlog_precomp(t, k % p);
    }
    dlog_precomp_clear(t);
}
void
dlog_bench(multilog_f flog, const ulong * p, const ulong * a, int np, int num)
{
    int i;
    ulong q;
    TIMEIT_ONCE_START
    for (i = 0; i < np; i ++)
        flog(p[i], a[i], num);
    TIMEIT_ONCE_STOP
}

int main()
{
    slong iter, k, nv, nref, r, nr;
    ulong minq, maxq;
    ulong * rand;
    int nbits, nl = 4;
    int vl[5] = { 1, 10, 100, 1000 , 5000};
    flint_rand_t state;

    flint_randinit(state);
    for (nbits = 10; nbits < 52; nbits += 10)
    {

        int i, np = 100;
        ulong p[100], a[100];
        for (i=0; i < np; i++)
        {
            p[i] = n_randprime(state, nbits, 0);
            a[i] = n_primitive_root_prime(p[i]);
        }

        for (i = 0; i < nl; i++)
        {    
            flint_printf("%wu * logs mod primes of size %wu.....\n", vl[i], nbits);

            if (nbits <= 20)
            { 
                flint_printf("table.......... ");
                fflush(stdout);
                dlog_bench(flog_table, p, a, np, vl[i]);
            }

            if (nbits <= 40 || vl[i] <= 10)
            {
                flint_printf("bsgs........... ");
                fflush(stdout);
                dlog_bench(flog_bsgs, p, a, np, vl[i]);
            }

            if (nbits <= 20 || vl[i] == 1)
            { 
                flint_printf("rho............ ");
                fflush(stdout);
                dlog_bench(flog_rho, p, a, np, vl[i]);
            }

            flint_printf("crt............ ");
            fflush(stdout);
            dlog_bench(flog_crt, p, a, np, vl[i]);

            flint_printf("generic........ ");
            fflush(stdout);
            dlog_bench(flog_gen, p, a, np, vl[i]);
        }
    }
    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
