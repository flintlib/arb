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

#include "acb_poly.h"
#include "bernoulli.h"

void _acb_poly_mullow_cpx(acb_ptr res, acb_srcptr src, long len, const acb_t c, long trunc, long prec);

void
_acb_poly_zeta_em_tail_naive(acb_ptr sum, const acb_t s, const acb_t Na, acb_srcptr Nasx, long M, long d, long prec)
{
    acb_ptr u, term;
    acb_t Na2, splus, rec;
    arb_t x;
    fmpz_t c;
    int aint;
    long r;

    BERNOULLI_ENSURE_CACHED(2 * M);

    u = _acb_vec_init(d);
    term = _acb_vec_init(d);
    acb_init(splus);
    acb_init(rec);
    acb_init(Na2);
    arb_init(x);
    fmpz_init(c);

    _acb_vec_zero(sum, d);

    /* u = 1/2 * Nasx */
    _acb_vec_scalar_mul_2exp_si(u, Nasx, d, -1L);

    /* term = u * (s+x) / (N+a) */
    _acb_poly_mullow_cpx(u, u, d, s, d, prec);
    _acb_vec_scalar_div(term, u, d, Na, prec);

    /* (N+a)^2 or 1/(N+a)^2 */
    acb_mul(Na2, Na, Na, prec);
    aint = acb_is_int(Na2);

    if (!aint)
        acb_inv(Na2, Na2, prec);

    for (r = 1; r <= M; r++)
    {
        /* printf("sum 2: %ld %ld\n", r, M); */

        /* sum += bernoulli number * term */
        arb_set_round_fmpz(x, fmpq_numref(bernoulli_cache + 2 * r), prec);
        arb_div_fmpz(x, x, fmpq_denref(bernoulli_cache + 2 * r), prec);

        _acb_vec_scalar_mul_arb(u, term, d, x, prec);
        _acb_vec_add(sum, sum, u, d, prec);

        /* multiply term by ((s+x)+2r-1)((s+x)+2r) / ((N+a)^2 * (2*r+1)*(2*r+2)) */
        acb_set(splus, s);
        arb_add_ui(acb_realref(splus), acb_realref(splus), 2*r-1, prec);
        _acb_poly_mullow_cpx(term, term, d, splus, d, prec);
        arb_add_ui(acb_realref(splus), acb_realref(splus), 1, prec);
        _acb_poly_mullow_cpx(term, term, d, splus, d, prec);

        /* TODO: combine with previous multiplication? */
        if (aint)
        {
            arb_mul_ui(x, acb_realref(Na2), 2*r+1, prec);
            arb_mul_ui(x, x, 2*r+2, prec);
            _acb_vec_scalar_div_arb(term, term, d, x, prec);
        }
        else
        {
            fmpz_set_ui(c, 2*r+1);
            fmpz_mul_ui(c, c, 2*r+2);
            acb_div_fmpz(rec, Na2, c, prec);
            _acb_vec_scalar_mul(term, term, d, rec, prec);
        }
    }

    _acb_vec_clear(u, d);
    _acb_vec_clear(term, d);
    acb_clear(splus);
    acb_clear(rec);
    acb_clear(Na2);
    arb_clear(x);
    fmpz_clear(c);
}

