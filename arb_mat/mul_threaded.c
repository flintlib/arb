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

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "arb_mat.h"
#include "pthread.h"

typedef struct
{
    arb_ptr * C;
    const arb_ptr * A;
    const arb_ptr * B;
    long ar0;
    long ar1;
    long bc0;
    long bc1;
    long br;
    long prec;
}
arb_mat_mul_arg_t;

void *
_arb_mat_mul_thread(void * arg_ptr)
{
    arb_mat_mul_arg_t arg = *((arb_mat_mul_arg_t *) arg_ptr);
    long i, j, k;

    for (i = arg.ar0; i < arg.ar1; i++)
    {
        for (j = arg.bc0; j < arg.bc1; j++)
        {
            arb_mul(arg.C[i] + j, arg.A[i] + 0, arg.B[0] + j, arg.prec);

            for (k = 1; k < arg.br; k++)
            {
                arb_addmul(arg.C[i] + j, arg.A[i] + k, arg.B[k] + j, arg.prec);
            }
        }
    }

    flint_cleanup();
    return NULL;
}

void
arb_mat_mul_threaded(arb_mat_t C, const arb_mat_t A, const arb_mat_t B, long prec)
{
    long ar, ac, br, bc, i, num_threads;
    pthread_t * threads;
    arb_mat_mul_arg_t * args;

    ar = arb_mat_nrows(A);
    ac = arb_mat_ncols(A);
    br = arb_mat_nrows(B);
    bc = arb_mat_ncols(B);

    if (ac != br || ar != arb_mat_nrows(C) || bc != arb_mat_ncols(C))
    {
        printf("arb_mat_mul_threaded: incompatible dimensions\n");
        abort();
    }

    if (br == 0)
    {
        arb_mat_zero(C);
        return;
    }

    if (A == C || B == C)
    {
        arb_mat_t T;
        arb_mat_init(T, ar, bc);
        arb_mat_mul_threaded(T, A, B, prec);
        arb_mat_swap(T, C);
        arb_mat_clear(T);
        return;
    }

    num_threads = flint_get_num_threads();
    threads = flint_malloc(sizeof(pthread_t) * num_threads);
    args = flint_malloc(sizeof(arb_mat_mul_arg_t) * num_threads);

    for (i = 0; i < num_threads; i++)
    {
        args[i].C = C->rows;
        args[i].A = A->rows;
        args[i].B = B->rows;

        if (ar >= bc)
        {
            args[i].ar0 = (ar * i) / num_threads;
            args[i].ar1 = (ar * (i + 1)) / num_threads;
            args[i].bc0 = 0;
            args[i].bc1 = bc;
        }
        else
        {
            args[i].ar0 = 0;
            args[i].ar1 = ar;
            args[i].bc0 = (bc * i) / num_threads;
            args[i].bc1 = (bc * (i + 1)) / num_threads;
        }

        args[i].br = br;
        args[i].prec = prec;
        pthread_create(&threads[i], NULL, _arb_mat_mul_thread, &args[i]);
    }

    for (i = 0; i < num_threads; i++)
    {
        pthread_join(threads[i], NULL);
    }

    flint_free(threads);
    flint_free(args);
}

