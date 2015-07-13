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

    Copyright (C) 2014 Fredrik Johansson

******************************************************************************/

#include "acb_hypgeom.h"

void
acb_hypgeom_pfq_sum(acb_t s, acb_t t, acb_srcptr a, long p,
    acb_srcptr b, long q, const acb_t z, long n, long prec)
{
    if (n > 4 && prec >= 128
        && _acb_vec_bits(a, p) * p + _acb_vec_bits(b, q) * q + 10 < prec / 2)
    {
        if (prec >= 256 && acb_bits(z) < prec * 0.01)
            acb_hypgeom_pfq_sum_bs(s, t, a, p, b, q, z, n, prec);
        else
            acb_hypgeom_pfq_sum_rs(s, t, a, p, b, q, z, n, prec);
    }
    else if (prec >= 1500 && n >= 30 + 100000 / (prec - 1000))
    {
        acb_hypgeom_pfq_sum_fme(s, t, a, p, b, q, z, n, prec);
    }
    else
    {
        acb_hypgeom_pfq_sum_forward(s, t, a, p, b, q, z, n, prec);
    }
}

void
acb_hypgeom_pfq_sum_invz(acb_t s, acb_t t, acb_srcptr a, long p,
    acb_srcptr b, long q, const acb_t z, const acb_t zinv, long n, long prec)
{
    if (n > 4 && prec >= 128
        && _acb_vec_bits(a, p) * p + _acb_vec_bits(b, q) * q + 10 < prec / 2)
    {
        if (prec >= 256 && acb_bits(zinv) < prec * 0.01)
            acb_hypgeom_pfq_sum_bs_invz(s, t, a, p, b, q, zinv, n, prec);
        else
            acb_hypgeom_pfq_sum_rs(s, t, a, p, b, q, z, n, prec);
    }
    else if (prec >= 1500 && n >= 30 + 100000 / (prec - 1000))
    {
        acb_hypgeom_pfq_sum_fme(s, t, a, p, b, q, z, n, prec);
    }
    else
    {
        acb_hypgeom_pfq_sum_forward(s, t, a, p, b, q, z, n, prec);
    }
}

