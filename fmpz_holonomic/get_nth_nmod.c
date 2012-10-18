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

mp_limb_t
fmpz_holonomic_get_nth_nmod(const fmpz_holonomic_t op, mp_srcptr initial, ulong n0, ulong n, nmod_t mod)
{
    long r = fmpz_holonomic_order(op);

    if (r == 0)
    {
        return 0;
    }
    else if (n < n0)
    {
        printf("not implemented\n");
        abort();
    }
    else if (n - n0 < r)
    {
        return initial[n - n0];
    }
    else
    {
        nmod_mat_t M;
        long i;
        mp_limb_t t, Q;
        nmod_mat_init(M, r, r, mod.n);

        fmpz_holonomic_forward_nmod_mat(M, &Q, op, n0, n - n0 - r + 1);

        t = 0;
        for (i = 0; i < r; i++)
        {
            NMOD_ADDMUL(t, M->rows[r - 1][i], initial[i], mod);
        }

        t = n_mulmod2_preinv(t, n_invmod(Q, mod.n), mod.n, mod.ninv);
        nmod_mat_clear(M);
        return t;
    }
}

