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

#include "fmprb_poly.h"

static __inline__ void fmpr_get_bot_exp(fmpz_t exp, const fmpr_t x)
{
    fmpz_set(exp, fmpr_expref(x));
}

static __inline__ void fmpr_get_top_exp(fmpz_t exp, const fmpr_t x)
{
    fmpz_add_ui(exp, fmpr_expref(x), fmpz_bits(fmpr_manref(x)));
}

void _fmpr_fmpz_vec_max_norm(fmpr_t norm, const fmpz * vec, long len, long prec)
{
    fmpr_set_fmpz(norm, vec + _fmpz_vec_height_index(vec, len));
    fmpr_set_round(norm, norm, prec, FMPR_RND_UP);
    fmpr_abs(norm, norm);
}

/* XXX: refactor this */
void fmprb_set_fmpz_2exp_round(fmprb_t y, const fmpz_t x, const fmpz_t exp, long prec)
{
    fmprb_set_fmpz(y, x);

    if (!fmpz_is_zero(x))
    {
        long r;

        fmpz_add(fmpr_expref(fmprb_midref(y)), fmpr_expref(fmprb_midref(y)), exp);

        r = fmpr_set_round(fmprb_midref(y), fmprb_midref(y), prec, FMPR_RND_DOWN);
        fmpr_set_error_result(fmprb_radref(y), fmprb_midref(y), r);
    }
}

int _fmprb_poly_mid_get_hull(fmpz_t bot_exp, fmpz_t top_exp, const fmprb_struct * A, long lenA)
{
    long i;
    fmpz_t t;
    int have_nonzero = 0;

    fmpz_init(t);
    fmpz_zero(bot_exp);
    fmpz_zero(top_exp);

    for (i = 0; i < lenA; i++)
    {
        if (fmpr_is_normal(fmprb_midref(A + i)))
        {
            if (!have_nonzero)
            {
                have_nonzero = 1;
                fmpr_get_bot_exp(bot_exp, fmprb_midref(A + i));
                fmpr_get_top_exp(top_exp, fmprb_midref(A + i));
            }
            else
            {
                fmpr_get_bot_exp(t, fmprb_midref(A + i));
                if (fmpz_cmp(t, bot_exp) < 0)
                    fmpz_swap(t, bot_exp);

                fmpr_get_top_exp(t, fmprb_midref(A + i));
                if (fmpz_cmp(t, top_exp) > 0)
                    fmpz_swap(t, top_exp);
            }
        }
        else if (!fmpr_is_zero(fmprb_midref(A + i)))
        {
            printf("exception: inf or nan encountered in polynomial\n");
            abort();
        }
    }

    fmpz_clear(t);
    return have_nonzero;
}

/* convert to an fmpz poly with a common exponent and coefficients
   at most prec bits, also bounding input error plus rounding error */
void _fmprb_poly_get_fmpz_poly_2exp(fmpr_t error, fmpz_t exp, fmpz  * coeffs,
                            const fmprb_struct * A, long lenA, long prec)
{
    fmpz_t top_exp, bot_exp;
    long shift;
    long i;
    int rounding;

    fmpz_init(top_exp);
    fmpz_init(bot_exp);

    if (!_fmprb_poly_mid_get_hull(bot_exp, top_exp, A, lenA))
    {
        fmpz_zero(exp);
        _fmpz_vec_zero(coeffs, lenA);
        fmpr_zero(error);
        for (i = 0; i < lenA; i++)
        {
            if (fmpr_cmp(fmprb_radref(A + i), error) > 0)
                fmpr_set(error, fmprb_radref(A + i));
        }

        return;   /* no need to clear fmpzs */
    }

    /* only take as much precision as necessary */
    shift = _fmpz_sub_small(top_exp, bot_exp);
    prec = FLINT_MIN(prec, shift);

    fmpz_sub_ui(exp, top_exp, prec);

    /* extract integer polynomial */
    rounding = 0;
    for (i = 0; i < lenA; i++)
    {
        if (fmpr_is_zero(fmprb_midref(A + i)))
        {
            fmpz_zero(coeffs + i);
        }
        else
        {
            shift = _fmpz_sub_small(fmpr_expref(fmprb_midref(A + i)), exp);

            if (shift >= 0)
            {
                fmpz_mul_2exp(coeffs + i,
                    fmpr_manref(fmprb_midref(A + i)), shift);
            }
            else
            {
                fmpz_tdiv_q_2exp(coeffs + i,
                    fmpr_manref(fmprb_midref(A + i)), -shift);
                rounding = 1;
            }
        }
    }

    fmpr_zero(error);

    /* compute maximum of input errors */
    for (i = 0; i < lenA; i++)
    {
        if (fmpr_cmp(fmprb_radref(A + i), error) > 0)
            fmpr_set(error, fmprb_radref(A + i));
    }

    /* add rounding error */
    if (rounding)
    {
        fmpr_t t;
        fmpr_init(t);

        fmpz_set_ui(fmpr_manref(t), 1UL);
        fmpz_set(fmpr_expref(t), exp);

        fmpr_add(error, error, t, FMPRB_RAD_PREC, FMPR_RND_UP);

        fmpr_clear(t);
    }

    fmpz_clear(top_exp);
}

void _fmprb_poly_mullow_ztrunc(fmprb_struct * C,
    const fmprb_struct * A, long lenA,
    const fmprb_struct * B, long lenB, long n, long prec)
{
    fmpz * Acoeffs, * Bcoeffs, * Ccoeffs;
    fmpz_t Aexp, Bexp, Cexp;
    fmpr_t Aerr, Berr, Anorm, Bnorm, err;
    long i;

    fmpz_init(Aexp);
    fmpz_init(Bexp);
    fmpz_init(Cexp);

    Acoeffs = _fmpz_vec_init(lenA);
    Bcoeffs = _fmpz_vec_init(lenB);
    Ccoeffs = _fmpz_vec_init(n);

    fmpr_init(Aerr);
    fmpr_init(Berr);
    fmpr_init(Anorm);
    fmpr_init(Bnorm);
    fmpr_init(err);

    _fmprb_poly_get_fmpz_poly_2exp(Aerr, Aexp, Acoeffs, A, lenA, prec);
    _fmprb_poly_get_fmpz_poly_2exp(Berr, Bexp, Bcoeffs, B, lenB, prec);

    /* main multiplication */
    if (lenA >= lenB)
        _fmpz_poly_mullow(Ccoeffs, Acoeffs, lenA, Bcoeffs, lenB, n);
    else
        _fmpz_poly_mullow(Ccoeffs, Bcoeffs, lenB, Acoeffs, lenA, n);

    fmpz_add(Cexp, Aexp, Bexp);

    /* cross-multiply error bounds: (A+r)(B+s) = AB + As + Br + rs */
    _fmpr_fmpz_vec_max_norm(Anorm, Acoeffs, lenA, FMPRB_RAD_PREC);
    fmpr_mul_2exp_fmpz(Anorm, Anorm, Aexp);
    _fmpr_fmpz_vec_max_norm(Bnorm, Bcoeffs, lenB, FMPRB_RAD_PREC);
    fmpr_mul_2exp_fmpz(Bnorm, Bnorm, Bexp);

    fmpr_mul(err, Aerr, Berr, FMPRB_RAD_PREC, FMPR_RND_UP);
    fmpr_addmul(err, Anorm, Berr, FMPRB_RAD_PREC, FMPR_RND_UP);
    fmpr_addmul(err, Bnorm, Aerr, FMPRB_RAD_PREC, FMPR_RND_UP);

    for (i = 0; i < n; i++)
    {
        fmprb_set_fmpz_2exp_round(C + i, Ccoeffs + i, Cexp, prec);

        /* there are at most (i+1) error terms for coefficient i */
        /* TODO: make this tight */
        fmpr_addmul_ui(fmprb_radref(C + i), err, i + 1,
            FMPRB_RAD_PREC, FMPR_RND_UP);
    }

    fmpr_clear(Aerr);
    fmpr_clear(Berr);
    fmpr_clear(Anorm);
    fmpr_clear(Bnorm);
    fmpr_clear(err);

    _fmpz_vec_clear(Acoeffs, lenA);
    _fmpz_vec_clear(Bcoeffs, lenB);
    _fmpz_vec_clear(Ccoeffs, n);

    fmpz_clear(Aexp);
    fmpz_clear(Bexp);
    fmpz_clear(Cexp);
}

void
fmprb_poly_mullow_ztrunc(fmprb_poly_t res, const fmprb_poly_t poly1,
              const fmprb_poly_t poly2, long n, long prec)
{
    long len_out;

    if (poly1->length == 0 || poly2->length == 0 || n == 0)
    {
        fmprb_poly_zero(res);
        return;
    }

    len_out = poly1->length + poly2->length - 1;
    if (n > len_out)
        n = len_out;

    fmprb_poly_fit_length(res, n);
    _fmprb_poly_mullow_ztrunc(res->coeffs, poly1->coeffs, poly1->length,
                                    poly2->coeffs, poly2->length, n, prec);
    _fmprb_poly_set_length(res, n);
    _fmprb_poly_normalise(res);
}
