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

void _fmprb_poly_get_fmpz_poly_2exp(fmpr_t error, fmpz_t exp, fmpz  * coeffs,
                            const fmprb_struct * A, long lenA, long prec);

void
_fmprb_vec_add_fmpz_vec_2exp(fmprb_struct * z,
    const fmprb_struct * x, const fmpz * y, const fmpz_t exp, long len, long prec)
{
    fmpr_t t;
    long i;

    fmpr_init(t);

    for (i = 0; i < len; i++)
    {
        fmpr_set_fmpz_2exp(t, y + i, exp);
        fmprb_add_fmpr(z + i, x + i, t, prec);
    }

    fmpr_clear(t);
}


/*
First a bunch of ugly code for bounding the errors.
This code should be refactored, and we should switch
to some kind of proper fmpr polynomial multiplication
when the lengths get extremely large

*/

#define ZERO      (LONG_MIN / 2)
#define MIN_SMALL (LONG_MIN / 8)
#define MAX_SMALL (-MIN_SMALL)

int read_exps(long * xmag, const fmpr_struct * x, long xlen)
{
    long i, mag;
    fmpz exp;

    for (i = 0; i < xlen; i++)
    {
        /* steps of two to read midpoints / radii */
        if (fmpr_is_zero(x + 2 * i))
        {
            xmag[i] = ZERO;
        }
        else
        {
            exp = *fmpr_expref(x + 2 * i);

            if (COEFF_IS_MPZ(exp))
                return 0;

            mag = exp + fmpz_bits(fmpr_manref(x + 2 * i));

            if (mag < MIN_SMALL || mag > MAX_SMALL)
                return 0;

            xmag[i] = mag;
        }
    }

    return 1;  /* all are small */
}

void
add_mulbound(fmprb_struct * z,
        const fmpr_struct * x, long xlen,
        const fmpr_struct * y, long ylen, long n)
{
    long i, j, k, a, b, count;
    long *xmag, *ymag;
    long mag, zmag;
    int xsmall, ysmall;
    fmpr_t err;

    xmag = flint_calloc(xlen + ylen, sizeof(long));
    ymag = xmag + xlen;

    fmpr_init(err);

    xsmall = read_exps(xmag, x, xlen);
    ysmall = read_exps(ymag, y, ylen);

    for (i = 0; i < n; i++)
    {
        count = 0;
        a = FLINT_MAX(0, i - ylen + 1);
        b = FLINT_MIN(i + 1, xlen);
        zmag = ZERO;
        count = b - a;

        if (xsmall && ysmall)
        {
            for (j = a; j < b; j++)
            {
                k = i - j;
                mag = xmag[j] + ymag[k];
                zmag = FLINT_MAX(zmag, mag);
            }

            if (zmag < MIN_SMALL + MIN_SMALL)
                count = 0;

            if (count != 0)
            {
                fmpr_set_ui(err, count);
                fmpr_mul_2exp_si(err, err, zmag);
                fmpr_add(fmprb_radref(z + i), fmprb_radref(z + i), err, FMPRB_RAD_PREC, FMPR_RND_UP);
            }
        }
        else
        {
            for (j = a; j < b; j++)
            {
                k = i - j;
                fmpr_mul(err, x + 2*j, y + 2*k, FMPRB_RAD_PREC, FMPR_RND_UP);
                fmpr_abs(err, err);
                fmpr_add(fmprb_radref(z + i), fmprb_radref(z + i), err, FMPRB_RAD_PREC, FMPR_RND_UP);
            }
        }
    }

    fmpr_clear(err);
    flint_free(xmag);
}


static int
is_exact(const fmprb_struct * x, long len)
{
    long i;

    for (i = 0; i < len; i++)
        if (!fmprb_is_exact(x + i))
            return 0;

    return 1;
}

static void
add_errors(fmprb_struct * z, const fmprb_struct * x, long xlen, const fmprb_struct * y, long ylen, long n)
{
    int xexact, yexact;

    xexact = is_exact(x, xlen);
    yexact = is_exact(y, ylen);

    if (!yexact)            add_mulbound(z, ((const fmpr_struct *) x), xlen, ((const fmpr_struct *) y) + 1, ylen, n);
    if (!xexact)            add_mulbound(z, ((const fmpr_struct *) x) + 1, xlen, ((const fmpr_struct *) y), ylen, n);
    if (!xexact && !yexact) add_mulbound(z, ((const fmpr_struct *) x) + 1, xlen, ((const fmpr_struct *) y) + 1, ylen, n);
}

/* largest block height = ALPHA*prec + BETA */
#define ALPHA 3.0
#define BETA 512
#define CLASSICAL_CUTOFF 6

static void
get_blocks(long * xblocks, const fmprb_struct * x, long xlen, long prec)
{
    fmpz_t top, bot, t, b, v, block_top, block_bot;
    long i, block;
    int in_zero;

    fmpz_init(top);
    fmpz_init(bot);
    fmpz_init(t);
    fmpz_init(b);
    fmpz_init(v);
    fmpz_init(block_top);
    fmpz_init(block_bot);

    xblocks[0] = 0;
    block = 0;
    in_zero = 1;

    for (i = 0; i < xlen; i++)
    {
        if (fmpr_is_special(fmprb_midref(x + i)))
            continue;

        fmpz_set(bot, fmpr_expref(fmprb_midref(x + i)));
        fmpz_add_ui(top, bot, fmpz_bits(fmpr_manref(fmprb_midref(x + i))));

        /* extend current block */
        if (in_zero)
        {
            fmpz_swap(block_top, top);
            fmpz_swap(block_bot, bot);
        }
        else
        {
            fmpz_max(t, top, block_top);
            fmpz_min(b, bot, block_bot);

            fmpz_sub(v, t, b);

            /* extend current block */
            if (prec == FMPR_PREC_EXACT || fmpz_cmp_ui(v, ALPHA * prec + BETA) < 0)
            {
                fmpz_swap(block_top, t);
                fmpz_swap(block_bot, b);
            }
            else  /* start new block */
            {
                block++;
                xblocks[block] = i;

                fmpz_swap(block_top, top);
                fmpz_swap(block_bot, bot);
            }
        }

        in_zero = 0;
    }

    xblocks[block + 1] = xlen;

    fmpz_clear(top);
    fmpz_clear(bot);
    fmpz_clear(t);
    fmpz_clear(b);
    fmpz_clear(v);
    fmpz_clear(block_top);
    fmpz_clear(block_bot);
}

static int
has_infnan(const fmprb_struct * x, long len)
{
    long i;

    for (i = 0; i < len; i++)
    {
        if (fmpr_is_nan(fmprb_midref(x + i)) || fmpr_is_inf(fmprb_midref(x + i))
            || fmpr_is_nan(fmprb_radref(x + i)) || fmpr_is_inf(fmprb_radref(x + i)))
        {
            return 1;
        }
    }

    return 0;
}

void
_fmprb_poly_mullow_block(fmprb_struct * z,
    const fmprb_struct * x, long xlen,
    const fmprb_struct * y, long ylen,
    long n, long prec)
{
    long *xblocks, *yblocks;
    long i, j, xp, yp, xl, yl, bn, r, s;
    fmprb_struct * tmp;
    fmpr_t t, error;
    fmpz_t xexp, yexp, zexp;
    fmpz *xz, *yz, *zz;

    if (has_infnan(x, xlen) || has_infnan(y, ylen))
    {
        _fmprb_poly_mullow_classical(z, x, xlen, y, ylen, n, prec);
        return;
    }

    xblocks = flint_malloc(sizeof(long) * (xlen + 1));
    yblocks = flint_malloc(sizeof(long) * (ylen + 1));

    fmpz_init(xexp);
    fmpz_init(yexp);
    fmpz_init(zexp);

    xz = _fmpz_vec_init(xlen);
    yz = _fmpz_vec_init(ylen);
    zz = _fmpz_vec_init(n);

    fmpr_init(error);
    fmpr_init(t);

    tmp = _fmprb_vec_init(n);

    /* split input polynomials into blocks with similar-size exponents */
    get_blocks(xblocks, x, xlen, prec);
    get_blocks(yblocks, y, ylen, prec);

    /* start with the zero polynomial plus the error */
    _fmprb_vec_zero(z, n);
    add_errors(z, x, xlen, y, ylen, n);

    for (i = 0; (xp = xblocks[i]) != xlen; i++)
    {
        for (j = 0; (yp = yblocks[j]) != ylen; j++)
        {
            if (xp + yp < n)
            {
                xl = xblocks[i + 1] - xp;
                yl = yblocks[j + 1] - yp;
                bn = FLINT_MIN(xl + yl - 1, n - xp - yp);
                xl = FLINT_MIN(xl, bn);
                yl = FLINT_MIN(yl, bn);

                /* classical multiplication; we just need the midpoints */
                if (bn < CLASSICAL_CUTOFF || xl < CLASSICAL_CUTOFF || yl < CLASSICAL_CUTOFF)
                {
                    for (r = 0; r < xl; r++)
                    {
                        for (s = 0; s < yl && xp + yp + r + s < n; s++)
                        {
                            long zi, ret;
                            zi = xp + yp + r + s;

                            fmpr_mul(t, fmprb_midref(x + xp + r), fmprb_midref(y + yp + s), FMPR_PREC_EXACT, FMPR_RND_DOWN);
                            ret = fmpr_add(fmprb_midref(z + zi), fmprb_midref(z + zi), t, prec, FMPR_RND_DOWN);

                            fmpr_add_error_result(fmprb_radref(z + zi), fmprb_radref(z + zi),
                                fmprb_midref(z + zi), ret, prec, FMPR_RND_UP);
                        }
                    }
                }
                else
                {
                    /* exact multiplication over Z[x]
                       todo: might want to verify that error==0 */

                    _fmprb_poly_get_fmpz_poly_2exp(error, xexp, xz, x + xp, xl, FMPR_PREC_EXACT);
                    _fmprb_poly_get_fmpz_poly_2exp(error, yexp, yz, y + yp, yl, FMPR_PREC_EXACT);
                    _fmpz_poly_mullow(zz, xz, xl, yz, yl, bn);
                    fmpz_add(zexp, xexp, yexp);
                    _fmprb_vec_add_fmpz_vec_2exp(z + xp + yp, z + xp + yp, zz, zexp, bn, prec);
                }
            }
        }
    }

    fmpz_clear(xexp);
    fmpz_clear(yexp);
    fmpz_clear(zexp);

    _fmpz_vec_clear(xz, xlen);
    _fmpz_vec_clear(yz, ylen);
    _fmpz_vec_clear(zz, n);

    fmpr_clear(error);
    fmpr_clear(t);

    _fmprb_vec_clear(tmp, n);
    flint_free(xblocks);
    flint_free(yblocks);
}

void
fmprb_poly_mullow_block(fmprb_poly_t res, const fmprb_poly_t poly1,
              const fmprb_poly_t poly2, long n, long prec)
{
    long xlen, ylen, zlen;

    xlen = poly1->length;
    ylen = poly2->length;

    if (xlen == 0 || ylen == 0 || n == 0)
    {
        fmprb_poly_zero(res);
        return;
    }

    xlen = FLINT_MIN(xlen, n);
    ylen = FLINT_MIN(ylen, n);
    zlen = FLINT_MIN(xlen + ylen - 1, n);

    if (res == poly1 || res == poly2)
    {
        fmprb_poly_t tmp;
        fmprb_poly_init2(tmp, zlen);
        _fmprb_poly_mullow_block(tmp->coeffs, poly1->coeffs, xlen,
            poly2->coeffs, ylen, zlen, prec);
        fmprb_poly_swap(res, tmp);
        fmprb_poly_clear(tmp);
    }
    else
    {
        fmprb_poly_fit_length(res, zlen);
        _fmprb_poly_mullow_block(res->coeffs, poly1->coeffs, xlen,
            poly2->coeffs, ylen, zlen, prec);
    }

    _fmprb_poly_set_length(res, zlen);
    _fmprb_poly_normalise(res);
}

