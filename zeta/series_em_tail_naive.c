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

    Copyright (C) 2012, 2013 Fredrik Johansson

******************************************************************************/

#include "zeta.h"
#include "fmpcb.h"
#include "fmpcb_poly.h"
#include "bernoulli.h"

void _fmpcb_poly_mullow_cpx(fmpcb_ptr res, fmpcb_srcptr src, long len, const fmpcb_t c, long trunc, long prec);


static __inline__ int
fmpcb_is_int(const fmpcb_t z)
{
    return fmprb_is_zero(fmpcb_imagref(z)) && fmprb_is_int(fmpcb_realref(z));
}

void
zeta_em_tail_naive(fmpcb_ptr sum, const fmpcb_t s, const fmpcb_t Na, fmpcb_srcptr Nasx, long M, long d, long prec)
{
    fmpcb_ptr u, term;
    fmpcb_t Na2, splus, rec;
    fmprb_t x;
    fmpz_t c;
    int aint;
    long r;

    BERNOULLI_ENSURE_CACHED(2 * M);

    u = _fmpcb_vec_init(d);
    term = _fmpcb_vec_init(d);
    fmpcb_init(splus);
    fmpcb_init(rec);
    fmpcb_init(Na2);
    fmprb_init(x);
    fmpz_init(c);

    _fmpcb_vec_zero(sum, d);

    /* u = 1/2 * Nasx */
    _fmpcb_vec_scalar_mul_2exp_si(u, Nasx, d, -1L);

    /* term = u * (s+x) / (N+a) */
    _fmpcb_poly_mullow_cpx(u, u, d, s, d, prec);
    _fmpcb_vec_scalar_div(term, u, d, Na, prec);

    /* (N+a)^2 or 1/(N+a)^2 */
    fmpcb_mul(Na2, Na, Na, prec);
    aint = fmpcb_is_int(Na2);

    if (!aint)
        fmpcb_inv(Na2, Na2, prec);

    for (r = 1; r <= M; r++)
    {
        /* printf("sum 2: %ld %ld\n", r, M); */

        /* sum += bernoulli number * term */
        fmprb_set_round_fmpz(x, fmpq_numref(bernoulli_cache + 2 * r), prec);
        fmprb_div_fmpz(x, x, fmpq_denref(bernoulli_cache + 2 * r), prec);

        _fmpcb_vec_scalar_mul_fmprb(u, term, d, x, prec);
        _fmpcb_vec_add(sum, sum, u, d, prec);

        /* multiply term by ((s+x)+2r-1)((s+x)+2r) / ((N+a)^2 * (2*r+1)*(2*r+2)) */
        fmpcb_set(splus, s);
        fmprb_add_ui(fmpcb_realref(splus), fmpcb_realref(splus), 2*r-1, prec);
        _fmpcb_poly_mullow_cpx(term, term, d, splus, d, prec);
        fmprb_add_ui(fmpcb_realref(splus), fmpcb_realref(splus), 1, prec);
        _fmpcb_poly_mullow_cpx(term, term, d, splus, d, prec);

        /* TODO: combine with previous multiplication? */
        if (aint)
        {
            fmprb_mul_ui(x, fmpcb_realref(Na2), 2*r+1, prec);
            fmprb_mul_ui(x, x, 2*r+2, prec);
            _fmpcb_vec_scalar_div_fmprb(term, term, d, x, prec);
        }
        else
        {
            fmpz_set_ui(c, 2*r+1);
            fmpz_mul_ui(c, c, 2*r+2);
            fmpcb_div_fmpz(rec, Na2, c, prec);
            _fmpcb_vec_scalar_mul(term, term, d, rec, prec);
        }
    }

    _fmpcb_vec_clear(u, d);
    _fmpcb_vec_clear(term, d);
    fmpcb_clear(splus);
    fmpcb_clear(rec);
    fmpcb_clear(Na2);
    fmprb_clear(x);
    fmpz_clear(c);
}

