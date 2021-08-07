/*
    Copyright (C) 2013-2014 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "partitions.h"
#include <pthread.h>

/* defined in flint*/
#define NUMBER_OF_SMALL_PARTITIONS 128
FLINT_DLL extern const unsigned int partitions_lookup[NUMBER_OF_SMALL_PARTITIONS];

slong partitions_hrr_needed_terms(double n);

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
hrr_sum_threaded(arb_t x, const fmpz_t n, slong N, int use_doubles)
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
partitions_fmpz_fmpz_hrr(fmpz_t p, const fmpz_t n, int use_doubles)
{
    arb_t x;
    arf_t bound;
    slong N;

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
        flint_printf("not unique!\n");
        arb_printd(x, 50);
        flint_printf("\n");
        flint_abort();
    }

    arb_clear(x);
    arf_clear(bound);
}

/* To compute p(n) mod 2^64. */
static void
partitions_vec(mp_ptr v, slong len)
{
    slong i, j, n;
    mp_limb_t p;

    for (n = 0; n < FLINT_MIN(len, NUMBER_OF_SMALL_PARTITIONS); n++)
        v[n] = partitions_lookup[n];

    for (n = NUMBER_OF_SMALL_PARTITIONS; n < len; n++)
    {
        p = 0;
        for (i = 1, j = 1; j <= n; j += 3 * i + 1, i++)
            p = v[n - j] - p;
        if ((i & 1) == 0)
            p = -p;
        for (i = 1, j = 2; j <= n; j += 3 * i + 2, i++)
            p = v[n - j] - p;
        if ((i & 1) != 0)
            p = -p;
        v[n] = p;
    }
}

/* The floor+vec method *requires* n <= 1498 for floor(p(n)/2^64)
   to be equal to floor(T/2^64). It is faster up to n ~= 1200.
   With doubles, it is faster up to n ~= 500. */
void
_partitions_fmpz_ui(fmpz_t res, ulong n, int use_doubles)
{
    if (n < NUMBER_OF_SMALL_PARTITIONS)
    {
        fmpz_set_ui(res, partitions_lookup[n]);
    }
    else if (FLINT_BITS == 64 && (n < 500 || (!use_doubles && n < 1200)))
    {
        mp_ptr tmp = flint_malloc((n + 1) * sizeof(mp_limb_t));

        if (n < 417)  /* p(n) < 2^64 */
        {
            partitions_vec(tmp, n + 1);
            fmpz_set_ui(res, tmp[n]);
        }
        else
        {
            arb_t x;
            arb_init(x);
            fmpz_set_ui(res, n);
            partitions_leading_fmpz(x, res, 4 * sqrt(n) - 50);
            arb_mul_2exp_si(x, x, -64);
            arb_floor(x, x, 4 * sqrt(n) - 50);

            if (arb_get_unique_fmpz(res, x))
            {
                fmpz_mul_2exp(res, res, 64);
                partitions_vec(tmp, n + 1);
                fmpz_add_ui(res, res, tmp[n]);
            }
            else
            {
                flint_printf("warning: failed at %wu\n", n);
                fmpz_set_ui(res, n);
                partitions_fmpz_fmpz_hrr(res, res, use_doubles);
            }
            arb_clear(x);
        }
        flint_free(tmp);
    }
    else
    {
        fmpz_set_ui(res, n);
        partitions_fmpz_fmpz_hrr(res, res, use_doubles);
    }
}

void
partitions_fmpz_fmpz(fmpz_t res, const fmpz_t n, int use_doubles)
{
    if (fmpz_cmp_ui(n, 2000) < 0)
    {
        if (fmpz_sgn(n) < 0)
            fmpz_zero(res);
        else
            _partitions_fmpz_ui(res, *n, use_doubles);
    }
    else
    {
        partitions_fmpz_fmpz_hrr(res, n, use_doubles);
    }
}

void
partitions_fmpz_ui(fmpz_t res, ulong n)
{
    _partitions_fmpz_ui(res, n, 0);
}

void
partitions_fmpz_ui_using_doubles(fmpz_t res, ulong n)
{
    _partitions_fmpz_ui(res, n, 1);
}

