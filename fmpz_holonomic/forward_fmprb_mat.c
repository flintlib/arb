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

#include "fmpz_holonomic.h"

void
_fmpz_holonomic_forward_bsplit_fmpz(fmpz_mat_t M, fmpz_t Q,
    const fmpz_holonomic_t op, long start, long a, long b);

void
_fmpz_holonomic_eval_companion_matrix_fmpz(fmpz_mat_t M, fmpz_t Q,
    const fmpz_holonomic_t op, long n);

void
_fmpz_holonomic_forward_bsplit_fmprb(fmprb_mat_t M, fmprb_t Q,
    const fmpz_holonomic_t op, long start, long a, long b, long prec)
{
    long r = fmpz_holonomic_order(op);

    if (b - a < 4)
    {
        fmpz_mat_t ZM;
        fmpz_t ZQ;

        fmpz_mat_init(ZM, r, r);
        fmpz_init(ZQ);

        _fmpz_holonomic_forward_bsplit_fmpz(ZM, ZQ, op, start, a, b);

        fmprb_mat_set_fmpz_mat(M, ZM);
        fmprb_set_fmpz(Q, ZQ);

        fmpz_mat_clear(ZM);
        fmpz_clear(ZQ);
    }
    else
    {
        fmprb_mat_t L, R;
        fmprb_t Q2;
        long m = a + (b - a) / 2;

        fmprb_mat_init(L, r, r);
        fmprb_mat_init(R, r, r);
        fmprb_init(Q2);

        _fmpz_holonomic_forward_bsplit_fmprb(L, Q, op, start, a, m, prec);
        _fmpz_holonomic_forward_bsplit_fmprb(R, Q2, op, start, m, b, prec);

        fmprb_mat_mul(M, R, L, prec);
        fmprb_mul(Q, Q, Q2, prec);

        fmprb_mat_clear(L);
        fmprb_mat_clear(R);
        fmprb_clear(Q2);
    }
}

void
fmpz_holonomic_forward_fmprb_mat(fmprb_mat_t M, fmprb_t Q,
    const fmpz_holonomic_t op, long start, long n, long prec)
{
    long r = fmpz_holonomic_order(op);

    if (r == 0 || n == 0)
    {
        fmprb_mat_one(M);
        fmprb_one(Q);
        return;
    }

    if (n < 0)
    {
        abort();
    }

    /* todo: hypergeom */

    if (fmpz_holonomic_seq_is_cfinite(op))
    {
        fmpz_mat_t ZM;
        fmpz_t ZQ;
        fmpz_mat_init(ZM, r, r);
        fmpz_init(ZQ);

        _fmpz_holonomic_eval_companion_matrix_fmpz(ZM, ZQ, op, 0);
        fmprb_mat_set_fmpz_mat(M, ZM);
        fmprb_set_fmpz(Q, ZQ);

        fmprb_mat_pow_ui(M, M, n, prec);
        fmprb_pow_ui(Q, Q, n, prec);

        fmpz_mat_clear(ZM);
        fmpz_clear(ZQ);
    }
    else
    {
        _fmpz_holonomic_forward_bsplit_fmprb(M, Q, op, start, 0, n, prec);
    }
}

