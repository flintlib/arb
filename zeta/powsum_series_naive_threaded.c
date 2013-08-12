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

#include <pthread.h>
#include "zeta.h"
#include "fmpcb.h"
#include "fmpcb_poly.h"

typedef struct
{
    fmpcb_ptr z;
    fmpcb_srcptr s;
    fmpcb_srcptr a;
    long n0;
    long n1;
    long len;
    long prec;
}
powsum_arg_t;

void _fmpcb_poly_fmpcb_invpow_cpx(fmpcb_ptr res,
    const fmpcb_t N, const fmpcb_t c, long trunc, long prec);

void *
_zeta_powsum_evaluator(void * arg_ptr)
{
    powsum_arg_t arg = *((powsum_arg_t *) arg_ptr);
    long k;

    fmpcb_ptr t;
    fmpcb_t c;

    t = _fmpcb_vec_init(arg.len);
    fmpcb_init(c);

    _fmpcb_vec_zero(arg.z, arg.len);

    for (k = arg.n0; k < arg.n1; k++)
    {
        fmpcb_add_ui(c, arg.a, k, arg.prec);
        _fmpcb_poly_fmpcb_invpow_cpx(t, c, arg.s, arg.len, arg.prec);
        _fmpcb_vec_add(arg.z, arg.z, t, arg.len, arg.prec);
    }

    _fmpcb_vec_clear(t, arg.len);
    fmpcb_clear(c);

    flint_cleanup();

    return NULL;
}

void
zeta_powsum_series_naive_threaded(fmpcb_ptr z,
    const fmpcb_t s, const fmpcb_t a, long n, long len, long prec)
{
    pthread_t * threads;
    powsum_arg_t * args;
    long i, num_threads;

    num_threads = flint_get_num_threads();

    threads = flint_malloc(sizeof(pthread_t) * num_threads);
    args = flint_malloc(sizeof(powsum_arg_t) * num_threads);

    for (i = 0; i < num_threads; i++)
    {
        args[i].z = _fmpcb_vec_init(len);
        args[i].s = s;
        args[i].a = a;
        args[i].n0 = (n * i) / num_threads;
        args[i].n1 = (n * (i + 1)) / num_threads;
        args[i].len = len;
        args[i].prec = prec;
        pthread_create(&threads[i], NULL, _zeta_powsum_evaluator, &args[i]);
    }

    for (i = 0; i < num_threads; i++)
    {
        pthread_join(threads[i], NULL);
    }

    _fmpcb_vec_zero(z, len);
    for (i = 0; i < num_threads; i++)
    {
        _fmpcb_vec_add(z, z, args[i].z, len, prec);
        _fmpcb_vec_clear(args[i].z, len);
    }

    flint_free(threads);
    flint_free(args);
}

