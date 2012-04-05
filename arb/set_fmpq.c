/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "arb.h"

void
_arb_set_fmpq(arb_t x, const fmpz_t p, const fmpz_t q)
{
    long pbits, qbits, shift;

    if (fmpz_is_one(q))
    {
        arb_set_fmpz(x, p);
        return;
    }

    pbits = fmpz_bits(p);
    qbits = fmpz_bits(q);
    shift = arb_prec(x) - (pbits - qbits);

    if (shift > 0)
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_mul_2exp(t, p, shift);
        fmpz_tdiv_q(arb_midref(x), t, q);
        fmpz_clear(t);
        fmpz_one(arb_radref(x));
        fmpz_set_si(arb_expref(x), -shift);
    }
    else
    {
        fmpz_tdiv_q(arb_midref(x), p, q);
        fmpz_zero(arb_expref(x));
        fmpz_one(arb_radref(x));
        _arb_normalise(x);
    }
}

void
arb_set_fmpq(arb_t x, const fmpq_t c)
{
    _arb_set_fmpq(x, fmpq_numref(c), fmpq_denref(c));
}
