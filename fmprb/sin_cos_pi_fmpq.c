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

#include "fmprb.h"

int
use_algebraic(const fmpz_t v, const fmpz_t w, long prec)
{
    fmpz q = *w;
    int r;

    if (COEFF_IS_MPZ(q))
        return 0;

    if (q <= 6)
        return 1;

    count_trailing_zeros(r, q);
    q >>= r;

    if (r >= 4 && prec < (r - 3) * 300)
        return 0;

    if (q > 1000)
        return 0;

    if (prec < 1500 + 150 * q)
        return 0;

    return 1;
}

void
_fmprb_sin_cos_pi_fmpq_oct(fmprb_t s, fmprb_t c,
            const fmpz_t v, const fmpz_t w, long prec)
{
    if (use_algebraic(v, w, prec))
    {
        _fmprb_sin_cos_pi_fmpq_algebraic(s, c, *v, *w, prec);
    }
    else
    {
        fmprb_const_pi(s, prec);
        fmprb_mul_fmpz(s, s, v, prec);
        fmprb_div_fmpz(s, s, w, prec);
        fmprb_sin_cos(s, c, s, prec);
    }
}

void
_fmprb_sin_pi_fmpq_oct(fmprb_t s, const fmpz_t v, const fmpz_t w, long prec)
{
    if (use_algebraic(v, w, prec))
    {
        _fmprb_sin_pi_fmpq_algebraic(s, *v, *w, prec);
    }
    else
    {
        fmprb_const_pi(s, prec);
        fmprb_mul_fmpz(s, s, v, prec);
        fmprb_div_fmpz(s, s, w, prec);
        fmprb_sin(s, s, prec);
    }
}

void
_fmprb_cos_pi_fmpq_oct(fmprb_t c, const fmpz_t v, const fmpz_t w, long prec)
{
    if (use_algebraic(v, w, prec))
    {
        _fmprb_cos_pi_fmpq_algebraic(c, *v, *w, prec);
    }
    else
    {
        fmprb_const_pi(c, prec);
        fmprb_mul_fmpz(c, c, v, prec);
        fmprb_div_fmpz(c, c, w, prec);
        fmprb_cos(c, c, prec);
    }
}

unsigned int
reduce_octant(fmpz_t v, fmpz_t w, const fmpq_t x)
{
    const fmpz * p = fmpq_numref(x);
    const fmpz * q = fmpq_denref(x);
    fmpz_t u;
    unsigned int octant;
    mp_bitcnt_t vval, wval;

    fmpz_init(u);

    fmpz_mul_2exp(w, p, 2);
    fmpz_fdiv_qr(w, v, w, q);
    octant = fmpz_fdiv_ui(w, 8);
    fmpz_mul_2exp(w, q, 2);

    if (octant % 2 != 0)
        fmpz_sub(v, q, v);

    vval = fmpz_val2(v);
    wval = fmpz_val2(w);
    vval = FLINT_MIN(vval, wval);

    if (vval != 0)
    {
        fmpz_tdiv_q_2exp(v, v, vval);
        fmpz_tdiv_q_2exp(w, w, vval);
    }

    fmpz_clear(u);

    return octant;
}

void
fmprb_sin_cos_pi_fmpq(fmprb_t s, fmprb_t c, const fmpq_t x, long prec)
{
    fmpz_t v, w;
    unsigned int octant;

    fmpz_init(v);
    fmpz_init(w);

    octant = reduce_octant(v, w, x);

    if ((octant + 1) % 4 < 2)
        _fmprb_sin_cos_pi_fmpq_oct(s, c, v, w, prec);
    else
        _fmprb_sin_cos_pi_fmpq_oct(c, s, v, w, prec);

    if ((octant + 6) % 8 < 4)  fmprb_neg(c, c);
    if (octant >= 4)           fmprb_neg(s, s);

    fmpz_clear(v);
    fmpz_clear(w);
}

void
fmprb_sin_pi_fmpq(fmprb_t s, const fmpq_t x, long prec)
{
    fmpz_t v, w;
    unsigned int octant;

    fmpz_init(v);
    fmpz_init(w);

    octant = reduce_octant(v, w, x);

    if ((octant + 1) % 4 < 2)
        _fmprb_sin_pi_fmpq_oct(s, v, w, prec);
    else
        _fmprb_cos_pi_fmpq_oct(s, v, w, prec);

    if (octant >= 4)
        fmprb_neg(s, s);

    fmpz_clear(v);
    fmpz_clear(w);
}

void
fmprb_cos_pi_fmpq(fmprb_t c, const fmpq_t x, long prec)
{
    fmpz_t v, w;
    unsigned int octant;

    fmpz_init(v);
    fmpz_init(w);

    octant = reduce_octant(v, w, x);

    if ((octant + 1) % 4 < 2)
        _fmprb_cos_pi_fmpq_oct(c, v, w, prec);
    else
        _fmprb_sin_pi_fmpq_oct(c, v, w, prec);

    if ((octant + 6) % 8 < 4)
        fmprb_neg(c, c);

    fmpz_clear(v);
    fmpz_clear(w);
}

