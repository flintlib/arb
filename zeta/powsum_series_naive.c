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

#include "zeta.h"
#include "fmpcb.h"
#include "fmpcb_poly.h"

void _fmpcb_poly_fmpcb_invpow_cpx(fmpcb_ptr res,
    const fmpcb_t N, const fmpcb_t c, long trunc, long prec);

void
zeta_powsum_series_naive(fmpcb_ptr z,
    const fmpcb_t s, const fmpcb_t a, long n, long len, long prec)
{
    long k;
    fmpcb_ptr t;
    fmpcb_t c;

    t = _fmpcb_vec_init(len);
    fmpcb_init(c);

    _fmpcb_vec_zero(z, len);

    for (k = 0; k < n; k++)
    {
        fmpcb_add_ui(c, a, k, prec);
        _fmpcb_poly_fmpcb_invpow_cpx(t, c, s, len, prec);
        _fmpcb_vec_add(z, z, t, len, prec);
    }

    _fmpcb_vec_clear(t, len);
    fmpcb_clear(c);
}

