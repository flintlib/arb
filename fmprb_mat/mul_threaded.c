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

#include "fmprb_mat.h"
#include "pthread.h"

typedef struct
{
    fmprb_ptr * C;
    const fmprb_ptr * A;
    const fmprb_ptr * B;
    long ar0;
    long ar1;
    long bc0;
    long bc1;
    long br;
    long prec;
}
fmprb_mat_mul_arg_t;

void *
_fmprb_mat_mul_thread(void * arg_ptr)
{
    fmprb_mat_mul_arg_t arg = *((fmprb_mat_mul_arg_t *) arg_ptr);
    long i, j, k;

    for (i = arg.ar0; i < arg.ar1; i++)
    {
        for (j = arg.bc0; j < arg.bc1; j++)
        {
            fmprb_mul(arg.C[i] + j, arg.A[i] + 0, arg.B[0] + j, arg.prec);

            for (k = 1; k < arg.br; k++)
            {
                fmprb_addmul(arg.C[i] + j, arg.A[i] + k, arg.B[k] + j, arg.prec);
            }
        }
    }

    flint_cleanup();
    return NULL;
}

void
fmprb_mat_mul_threaded(fmprb_mat_t C, const fmprb_mat_t A, const fmprb_mat_t B, long prec)
{
    long ar, ac, br, bc, i, num_threads;
    pthread_t * threads;
    fmprb_mat_mul_arg_t * args;

    ar = fmprb_mat_nrows(A);
    ac = fmprb_mat_ncols(A);
    br = fmprb_mat_nrows(B);
    bc = fmprb_mat_ncols(B);

    if (ac != br || ar != fmprb_mat_nrows(C) || bc != fmprb_mat_ncols(C))
    {
        printf("fmprb_mat_mul_threaded: incompatible dimensions\n");
        abort();
    }

    if (br == 0)
    {
        fmprb_mat_zero(C);
        return;
    }

    if (A == C || B == C)
    {
        fmprb_mat_t T;
        fmprb_mat_init(T, ar, bc);
        fmprb_mat_mul_threaded(T, A, B, prec);
        fmprb_mat_swap(T, C);
        fmprb_mat_clear(T);
        return;
    }

    num_threads = flint_get_num_threads();
    threads = flint_malloc(sizeof(pthread_t) * num_threads);
    args = flint_malloc(sizeof(fmprb_mat_mul_arg_t) * num_threads);

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
        pthread_create(&threads[i], NULL, _fmprb_mat_mul_thread, &args[i]);
    }

    for (i = 0; i < num_threads; i++)
    {
        pthread_join(threads[i], NULL);
    }

    flint_free(threads);
    flint_free(args);
}

