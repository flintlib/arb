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

    Copyright (C) 2013-2014 Fredrik Johansson

******************************************************************************/

#include <pthread.h>
#include "partitions.h"

/* defined in flint*/
#define NUMBER_OF_SMALL_PARTITIONS 128
extern const unsigned int partitions_lookup[NUMBER_OF_SMALL_PARTITIONS];

long partitions_hrr_needed_terms(double n);

typedef struct
{
    arb_ptr x;
    fmpz * n;
    ulong N0;
    ulong N;
    int use_doubles;
}
worker_arg_t;

static void *
worker(void * arg_ptr)
{
    worker_arg_t arg = *((worker_arg_t *) arg_ptr);
    partitions_hrr_sum_arb(arg.x, arg.n, arg.N0, arg.N, arg.use_doubles);
    flint_cleanup();
    return NULL;
}

/* TODO: set number of threads in child threads, for future
   multithreaded evaluation of single terms */
static void
hrr_sum_threaded(arb_t x, const fmpz_t n, long N, int use_doubles)
{
    arb_t y;
    pthread_t threads[2];
    worker_arg_t args[2];

    arb_init(y);

    args[0].x = x;
    args[0].n = (fmpz *) n;
    args[0].N0 = 1;
    args[0].N = 16;
    args[0].use_doubles = use_doubles;

    args[1].x = y;
    args[1].n = (fmpz *) n;
    args[1].N0 = 17;
    args[1].N = N;
    args[1].use_doubles = use_doubles;

    pthread_create(&threads[0], NULL, worker, &args[0]);
    pthread_create(&threads[1], NULL, worker, &args[1]);

    pthread_join(threads[0], NULL);
    pthread_join(threads[1], NULL);

    arb_add(x, x, y, ARF_PREC_EXACT);

    arb_clear(y);
}

void
partitions_fmpz_fmpz(fmpz_t p, const fmpz_t n, int use_doubles)
{
    if (fmpz_cmp_ui(n, NUMBER_OF_SMALL_PARTITIONS) < 0)
    {
        if (fmpz_sgn(n) < 0)
            fmpz_zero(p);
        else
            fmpz_set_ui(p, partitions_lookup[fmpz_get_ui(n)]);
    }
    else
    {
        arb_t x;
        arf_t bound;
        long N;

        arb_init(x);
        arf_init(bound);

        N = partitions_hrr_needed_terms(fmpz_get_d(n));

        if (fmpz_cmp_ui(n, 4e8) >= 0 && flint_get_num_threads() > 1)
        {
            hrr_sum_threaded(x, n, N, use_doubles);
        }
        else
        {
            partitions_hrr_sum_arb(x, n, 1, N, use_doubles);
        }

        partitions_rademacher_bound(bound, n, N);
        arb_add_error_arf(x, bound);

        if (!arb_get_unique_fmpz(p, x))
        {
            printf("not unique!\n");
            arb_printd(x, 50);
            printf("\n");
            abort();
        }

        arb_clear(x);
        arf_clear(bound);
    }
}

void
partitions_fmpz_ui(fmpz_t p, ulong n)
{
    if (n < NUMBER_OF_SMALL_PARTITIONS)
    {
        fmpz_set_ui(p, partitions_lookup[n]);
    }
    else
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_set_ui(t, n);
        partitions_fmpz_fmpz(p, t, 0);
        fmpz_clear(t);
    }
}

void
partitions_fmpz_ui_using_doubles(fmpz_t p, ulong n)
{
    if (n < NUMBER_OF_SMALL_PARTITIONS)
    {
        fmpz_set_ui(p, partitions_lookup[n]);
    }
    else
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_set_ui(t, n);
        partitions_fmpz_fmpz(p, t, 1);
        fmpz_clear(t);
    }
}

