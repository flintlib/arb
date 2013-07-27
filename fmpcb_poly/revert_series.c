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

#include "fmpcb_poly.h"

void
_fmpcb_poly_revert_series(fmpcb_ptr Qinv,
    fmpcb_srcptr Q, long n, long prec)
{
    _fmpcb_poly_revert_series_lagrange_fast(Qinv, Q, n, prec);
}

void
fmpcb_poly_revert_series(fmpcb_poly_t Qinv,
                                    const fmpcb_poly_t Q, long n, long prec)
{
    fmpcb_ptr Qcopy;
    int Qalloc;
    long Qlen = Q->length;

    if (Q->length < 2 || !fmpcb_is_zero(Q->coeffs)
                      || fmpcb_contains_zero(Q->coeffs + 1))
    {
        printf("Exception (fmpcb_poly_revert_series). Input must \n"
               "have zero constant term and nonzero coefficient of x^1.\n");
        abort();
    }

    if (n < 2)
    {
        fmpcb_poly_zero(Qinv);
        return;
    }

    if (Qlen >= n)
    {
        Qcopy = Q->coeffs;
        Qalloc = 0;
    }
    else
    {
        long i;
        Qcopy = _fmpcb_vec_init(n);
        for (i = 0; i < Qlen; i++)
            Qcopy[i] = Q->coeffs[i];
        Qalloc = 1;
    }

    if (Qinv != Q)
    {
        fmpcb_poly_fit_length(Qinv, n);
        _fmpcb_poly_revert_series(Qinv->coeffs, Qcopy, n, prec);
    }
    else
    {
        fmpcb_poly_t t;
        fmpcb_poly_init2(t, n);
        _fmpcb_poly_revert_series(t->coeffs, Qcopy, n, prec);
        fmpcb_poly_swap(Qinv, t);
        fmpcb_poly_clear(t);
    }

    _fmpcb_poly_set_length(Qinv, n);
    _fmpcb_poly_normalise(Qinv);

    if (Qalloc)
        flint_free(Qcopy);
}

