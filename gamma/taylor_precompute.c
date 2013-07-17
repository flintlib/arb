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

#include "gamma.h"
#include "fmprb_poly.h"

TLS_PREFIX fmprb_struct * gamma_taylor_coeffs = NULL;

TLS_PREFIX long gamma_taylor_prec = 0;

TLS_PREFIX long gamma_taylor_num = 0;

void
gamma_taylor_precompute(long num, long prec)
{
    long wp;

    if (gamma_taylor_num < num || gamma_taylor_prec < prec)
    {
        fmprb_poly_t A;
        long i;

        _fmprb_vec_clear(gamma_taylor_coeffs, gamma_taylor_num);

        num = FLINT_MAX(gamma_taylor_num * 1.5, num);
        prec = FLINT_MAX(gamma_taylor_prec * 1.5, prec);
        prec = FLINT_MAX(prec, 64);

        wp = prec * 1.01 + 32;

        /* TODO: cleanup */
        fmprb_poly_init(A);
        fmprb_poly_log_gamma_series(A, num, wp);
        fmprb_poly_neg(A, A);
        fmprb_poly_exp_series(A, A, num, wp);

        for (i = 1; i < num; i += 2)
            fmprb_neg(A->coeffs + i, A->coeffs + i);

        gamma_taylor_coeffs = A->coeffs;
        gamma_taylor_num = A->length;
        gamma_taylor_prec = prec;
    }
}

