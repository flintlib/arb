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

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "bernoulli.h"

long bernoulli_cache_num = 0;

fmpq * bernoulli_cache = NULL;

void
bernoulli_cache_compute(long n)
{
    if (bernoulli_cache_num < n)
    {
        long i, new_num;
        bernoulli_rev_t iter;

        new_num = FLINT_MAX(bernoulli_cache_num * 2, n);
        new_num = n;

        bernoulli_cache = flint_realloc(bernoulli_cache, new_num * sizeof(fmpq));
        for (i = bernoulli_cache_num; i < new_num; i++)
            fmpq_init(bernoulli_cache + i);

        i = new_num - 1;
        i -= (i % 2);
        bernoulli_rev_init(iter, i);
        for ( ; i >= 0; i -= 2)
            bernoulli_rev_next(fmpq_numref(bernoulli_cache + i),
                fmpq_denref(bernoulli_cache + i), iter);
        bernoulli_rev_clear(iter);

        if (new_num > 1)
            fmpq_set_si(bernoulli_cache + 1, -1, 2);

        bernoulli_cache_num = new_num;
    }
}

