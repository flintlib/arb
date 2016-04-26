/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "bernoulli.h"

TLS_PREFIX slong bernoulli_cache_num = 0;

TLS_PREFIX fmpq * bernoulli_cache = NULL;

void
bernoulli_cleanup(void)
{
    slong i;

    for (i = 0; i < bernoulli_cache_num; i++)
        fmpq_clear(bernoulli_cache + i);

    flint_free(bernoulli_cache);
    bernoulli_cache = NULL;
    bernoulli_cache_num = 0;
}

void
bernoulli_cache_compute(slong n)
{
    if (bernoulli_cache_num < n)
    {
        slong i, new_num;
        bernoulli_rev_t iter;

        if (bernoulli_cache_num == 0)
        {
            flint_register_cleanup_function(bernoulli_cleanup);
        }

        new_num = FLINT_MAX(bernoulli_cache_num + 128, n);

        bernoulli_cache = flint_realloc(bernoulli_cache, new_num * sizeof(fmpq));
        for (i = bernoulli_cache_num; i < new_num; i++)
            fmpq_init(bernoulli_cache + i);

        i = new_num - 1;
        i -= (i % 2);
        bernoulli_rev_init(iter, i);
        for ( ; i >= bernoulli_cache_num; i -= 2)
        {
            bernoulli_rev_next(fmpq_numref(bernoulli_cache + i),
                fmpq_denref(bernoulli_cache + i), iter);
        }
        bernoulli_rev_clear(iter);

        if (new_num > 1)
            fmpq_set_si(bernoulli_cache + 1, -1, 2);

        bernoulli_cache_num = new_num;
    }
}

