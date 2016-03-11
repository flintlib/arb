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

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("dlog mod p....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        dlog_table_t table;
        dlog_bsgs_t bsgs;
        dlog_rho_t rho;
        dlog_precomp_t pre1, pre100;
        ulong p, a, k;

        p = n_randprime(state, 15, 0);
        a = n_primitive_root_prime(p);

        dlog_table_init(table, a, p);
        dlog_bsgs_init(bsgs, a, p, p-1, ceil(sqrt((double)p)));
        dlog_rho_init(rho, a, p, p-1);
        dlog_precomp_n_init(pre1, a, p, p-1, 1);
        dlog_precomp_n_init(pre100, a, p, p-1, 100);

        for (k = 1; k < 100 && k < p; k++)
        {
            ulong l1, l2, l3, l4, l5;
            l1 = dlog_table(table, k);
            l2 = dlog_bsgs(bsgs, k);
            l3 = dlog_rho(rho, k);
            l4 = dlog_precomp(pre1, k);
            l5 = dlog_precomp(pre100, k);
            if (l1 != l2 || l1 != l3 || l1 != l4 || l1 != l5)
            {
                flint_printf("FAIL: log(%wu,%wu) mod %wu: [%wu, %wu, %wu, %wu, %wu]\n",
                        k, a, p, l1, l2, l3, l4, l5);
                abort();
            }
        }
        dlog_table_clear(table);
        dlog_bsgs_clear(bsgs);
        dlog_rho_clear(rho);
        dlog_precomp_clear(pre1);
        dlog_precomp_clear(pre100);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
