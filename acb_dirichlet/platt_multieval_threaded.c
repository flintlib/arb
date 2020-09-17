/*
    Copyright (C) 2020 Rudolph
    Copyright (C) 2020 D.H.J. Polymath

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"
#include "pthread.h"

void _platt_smk(acb_ptr table, const arb_t t0, slong A, slong B,
        slong Jstart, slong Jstop, slong K, slong prec);

void _acb_dirichlet_platt_multieval(arb_ptr out, acb_srcptr S_table,
        const arb_t t0, slong A, slong B, const arb_t h, slong J,
        slong K, slong sigma, slong prec);

typedef struct
{
    acb_ptr threadtable;
    arb_srcptr t0;
    slong A;
    slong B;
    slong K;
    slong start;
    slong stop;
    slong prec;
}
platt_smk_arg_t;

static void *
_platt_smk_thread(void * arg_ptr)
{
    platt_smk_arg_t *p = (platt_smk_arg_t *) arg_ptr;
    _platt_smk(p->threadtable, p->t0, p->A, p->B,
               p->start, p->stop, p->K, p->prec);
    flint_cleanup();
    return NULL;
}

void
acb_dirichlet_platt_multieval_threaded(arb_ptr out, const fmpz_t T, slong A,
        slong B, const arb_t h, slong J, slong K, slong sigma, slong prec)
{
    slong i, num_threads, N, threadtasks;
    pthread_t * threads;
    platt_smk_arg_t * args;
    acb_ptr S;
    arb_t t0;

    N = A*B;
    num_threads = flint_get_num_threads();
    threads = flint_malloc(sizeof(pthread_t) * num_threads);
    args = flint_malloc(sizeof(platt_smk_arg_t) * num_threads);
    threadtasks = (J+num_threads-1)/num_threads;

    arb_init(t0);
    arb_set_fmpz(t0, T);

    S =  _acb_vec_init(K*N);
    for (i = 0; i < num_threads; i++)
    {
        args[i].threadtable = _acb_vec_init(K*N);
        args[i].t0 = t0;
        args[i].A = A;
        args[i].B = B;
        args[i].K = K;
        args[i].start = i*threadtasks + 1;
        args[i].stop = (i+1)*threadtasks;;
        args[i].prec = prec;
    }
    args[num_threads-1].stop = J;

    for (i = 0; i < num_threads; i++)
    {
        pthread_create(&threads[i], NULL, _platt_smk_thread, &args[i]);
    }

    for (i = 0; i < num_threads; i++)
    {
        pthread_join(threads[i], NULL);
    }

    _acb_vec_zero(S, K*N);
    for (i = 0; i < num_threads; i++)
    {
        _acb_vec_add(S, S, args[i].threadtable, K*N, prec);
        _acb_vec_clear(args[i].threadtable, K*N);
    }

    _acb_dirichlet_platt_multieval(out, S, t0, A, B, h, J, K, sigma, prec);

    arb_clear(t0);
    _acb_vec_clear(S, K*N);
}
