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

    flint_printf("modpe....");
    fflush(stdout);
    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        ulong p, e, pe, a;

        if (iter == 0)
        {
            p = 2;
            pe = 8;
            a = 5;
            e = 3;
        }
        else
        {
            p = pe = n_randprime(state, 10, 0);
            a = (p == 40487) ? 10 : n_primitive_root_prime(p);
            e = 1;
        }

        for (; pe < UWORD(1000000000000); pe *= p, e++)
        {
            ulong k, phi;
            nmod_t mod;

            dlog_modpe_t modpe;

            nmod_init(&mod, pe);
            phi = (p == 2) ? pe / 4 : pe - pe / p;

            dlog_modpe_init(modpe, a, p, e, pe, 10);

            for (k = 0; k < 100 && k < p; k++)
            {
                ulong l, b, x;

                l = n_randint(state, phi);
                b = nmod_pow_ui(a, l, mod);

                if ((x = dlog_modpe(modpe, b)) != l)
                {
                    flint_printf("FAIL modpe: %wu^%wu = %wu [%wu^%wu]\n\n",
                            a, l, b, p, e);
                    flint_printf("modpe returned %wu\n\n", x);
                    abort();
                }
            }

            dlog_modpe_clear(modpe);
        }
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
