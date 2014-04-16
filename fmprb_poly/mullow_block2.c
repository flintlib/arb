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

#include <math.h>
#include "fmprb_poly.h"

void _fmprb_poly_get_scale(fmpz_t scale, fmprb_srcptr x, long xlen,
                                         fmprb_srcptr y, long ylen);

static int
_fmprb_vec_is_finite(fmprb_srcptr x, long len)
{
    long i;

    for (i = 0; i < len; i++)
        if (!fmprb_is_finite(x + i))
            return 0;

    return 1;
}

/*static __inline__ */ void
fmpr_add_abs_ubound(fmpr_t z, const fmpr_t x, const fmpr_t y, long prec)
{
    if (fmpr_sgn(x) >= 0)
    {
        if (fmpr_sgn(y) >= 0)
            fmpr_add(z, x, y, prec, FMPR_RND_UP);
        else
            fmpr_sub(z, x, y, prec, FMPR_RND_UP);
    }
    else
    {
        if (fmpr_sgn(y) >= 0)
            fmpr_sub(z, y, x, prec, FMPR_RND_UP);
        else
        {
            fmpr_add(z, x, y, prec, FMPR_RND_UP);
            fmpr_neg(z, z);
        }
    }
}

static __inline__ void
fmpr_abs_round(fmpr_t z, const fmpr_t x, long prec, fmpr_rnd_t rnd)
{
    if (fmpr_sgn(x) >= 0)
        fmpr_set_round(z, x, prec, rnd);
    else
        fmpr_neg_round(z, x, prec, rnd);
}

/* Break fmpr vector into same-exponent blocks where the largest block
   has a height of at most ALPHA*prec + BETA bits. These are just
   tuning parameters. Note that ALPHA * FMPRB_RAD_PREC + BETA
   should be smaller than DOUBLE_BLOCK_MAX_HEIGHT if we want to use
   doubles for error bounding. */
#define ALPHA 3.0
#define BETA 512


/* Maximum length of block for which we use double multiplication
   (for longer blocks, we use fmpz_poly multiplication). This is essentially
   just a tuning parameter, but note that it must be considered when
   compensating for rounding error below. */
#define DOUBLE_BLOCK_MAX_LENGTH 1000

/* Computing a dot product of length DOUBLE_BLOCK_MAX_LENGTH involving
   only nonnegative numbers, and then multiplying by this factor, must give
   an upper bound for the exact dot product (we can assume that no
   overflow or underflow occurs). The following is certainly
   sufficient, but it would be nice to include a formal proof here. */
#define DOUBLE_ROUNDING_FACTOR (1.0 + 1e-9)

/* Maximum height for which we use double multiplication. Since the dynamic
   exponent range of doubles is about +/- 1024, this must be less than about
   1024 (to allow the product of two numbers). This must also
   account for adding FMPRB_RAD_PREC bits. */
#define DOUBLE_BLOCK_MAX_HEIGHT 800

/* We divide coefficients by 2^DOUBLE_BLOCK_SHIFT when converting them to
   doubles, in order to use the whole exponent range. Note that this means
   numbers of size (2^(-DOUBLE_BLOCK_SHIFT))^2 must not underflow. */
#define DOUBLE_BLOCK_SHIFT (DOUBLE_BLOCK_MAX_HEIGHT / 2)


/* Converts fmpr vector to a vector of blocks, where each block is
   a vector of fmpz integers on a common exponent.
   Optionally, we also generate doubles on the same common exponent
   (minus DOUBLE_BLOCK_SHIFT), where this can be done exactly. */

static __inline__ int           /* returns new can_use_doubles status */
_fmpr_vec_get_fmpz_2exp_blocks
(
    fmpz * coeffs,            /* output fmpz coefficients */
    double * dblcoeffs,       /* output double coefficients (optional) */
    fmpz * exps,              /* common exponent of each block */
    long * blocks,            /* start positions of blocks (plus end marker) */
    const fmpz_t scale,       /* compose input poly by x -> x/2^scale */
    fmpr_srcptr x,            /* first in vector of input coefficient */
    long len,                 /* number of input coefficients */
    long step,                /* step length to read input coefficients*/
    long prec,                /* prec */
    int can_use_doubles
)
{
    fmpz_t top, bot, t, b, v, block_top, block_bot;
    long i, j, s, block, bits, maxheight;
    int in_zero;

    fmpz_init(top);
    fmpz_init(bot);
    fmpz_init(t);
    fmpz_init(b);
    fmpz_init(v);
    fmpz_init(block_top);
    fmpz_init(block_bot);

    blocks[0] = 0;
    block = 0;
    in_zero = 1;

    if (prec == FMPR_PREC_EXACT)
        maxheight = FMPR_PREC_EXACT;
    else
        maxheight = ALPHA * prec + BETA;

    can_use_doubles = can_use_doubles && (maxheight <= DOUBLE_BLOCK_MAX_HEIGHT);

    for (i = 0; i < len; i++)
    {
        /* Skip (must be zero, since we assume there are no Infs/NaNs). */
        if (fmpr_is_special(x + i * step))
            continue;

        /* Bottom and top exponent of current number */
        fmpz_set(bot, fmpr_expref(x + i * step));
        fmpz_submul_ui(bot, scale, i);
        bits = fmpz_bits(fmpr_manref(x + i * step));
        fmpz_add_ui(top, bot, bits);

        can_use_doubles = can_use_doubles && (bits <= FMPRB_RAD_PREC);

        /* Extend current block. */
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
            if (fmpz_cmp_ui(v, maxheight) < 0)
            {
                fmpz_swap(block_top, t);
                fmpz_swap(block_bot, b);
            }
            else  /* start new block */
            {
                /* write exponent for previous block */
                fmpz_set(exps + block, block_bot);

                block++;
                blocks[block] = i;

                fmpz_swap(block_top, top);
                fmpz_swap(block_bot, bot);
            }
        }

        in_zero = 0;
    }

    /* write exponent for last block */
    fmpz_set(exps + block, block_bot);

    /* end marker */
    blocks[block + 1] = len;

    /* write the block data */
    for (i = 0; blocks[i] != len; i++)
    {
        for (j = blocks[i]; j < blocks[i + 1]; j++)
        {
            if (fmpr_is_special(x + j * step))
            {
                fmpz_zero(coeffs + j);

                if (can_use_doubles)
                    dblcoeffs[j] = 0.0;
            }
            else
            {
                /* Divide coefficient by 2^(scale * j) */
                fmpz_mul_ui(t, scale, j);
                fmpz_sub(t, fmpr_expref(x + j * step), t);
                s = _fmpz_sub_small(t, exps + i);

                if (s < 0) abort(); /* Bug catcher */

                fmpz_mul_2exp(coeffs + j, fmpr_manref(x + j * step), s);

                if (can_use_doubles)
                {
                    double c = *fmpr_manref(x + j * step);
                    c = ldexp(c, s - DOUBLE_BLOCK_SHIFT);

                    if (c < 1e-150 || c > 1e150) /* Bug catcher */
                        abort();

                    dblcoeffs[j] = c;
                }
            }
        }
    }

    fmpz_clear(top);
    fmpz_clear(bot);
    fmpz_clear(t);
    fmpz_clear(b);
    fmpz_clear(v);
    fmpz_clear(block_top);
    fmpz_clear(block_bot);

    return can_use_doubles;
}

static __inline__ void
_fmprb_poly_addmullow_rad(fmprb_ptr z, fmpz * zz,
    const fmpz * xz, const double * xdbl, const fmpz * xexps,
    const long * xblocks, long xlen,
    const fmpz * yz, const double * ydbl, const fmpz * yexps,
    const long * yblocks, long ylen,
    long n, int can_use_doubles)
{
    long i, j, k, ii, xp, yp, xl, yl, bn;
    fmpz_t zexp;
    fmpr_t t;

    fmpz_init(zexp);
    fmpr_init(t);

    for (i = 0; (xp = xblocks[i]) != xlen; i++)
    {
        for (j = 0; (yp = yblocks[j]) != ylen; j++)
        {
            if (xp + yp >= n)
                continue;

            xl = xblocks[i + 1] - xp;
            yl = yblocks[j + 1] - yp;
            bn = FLINT_MIN(xl + yl - 1, n - xp - yp);
            xl = FLINT_MIN(xl, bn);
            yl = FLINT_MIN(yl, bn);

            fmpz_add_inline(zexp, xexps + i, yexps + j);

            if (can_use_doubles && xl > 1 && yl > 1 &&
                (xl < DOUBLE_BLOCK_MAX_LENGTH || yl < DOUBLE_BLOCK_MAX_LENGTH))
            {
                fmpz_add_ui(zexp, zexp, 2 * DOUBLE_BLOCK_SHIFT);

                for (k = 0; k < bn; k++)
                {
                    /* Classical multiplication (may round down!) */
                    double ss = 0.0;

                    for (ii = FLINT_MAX(0, k - yl + 1);
                        ii <= FLINT_MIN(xl - 1, k); ii++)
                    {
                        ss += xdbl[xp + ii] * ydbl[yp + k - ii];
                    }

                    /* Compensate for rounding error */
                    ss *= DOUBLE_ROUNDING_FACTOR;

                    fmpr_set_d(t, ss);
                    fmpr_mul_2exp_fmpz(t, t, zexp);
                    fmpr_add(fmprb_radref(z + xp + yp + k),
                        fmprb_radref(z + xp + yp + k), t,
                        FMPRB_RAD_PREC, FMPR_RND_UP);
                }
            }
            else
            {
                if (xl >= yl)
                    _fmpz_poly_mullow(zz, xz + xp, xl, yz + yp, yl, bn);
                else
                    _fmpz_poly_mullow(zz, yz + yp, yl, xz + xp, xl, bn);

                for (k = 0; k < bn; k++)
                {
                    fmpr_set_round_fmpz_2exp(t, zz + k, zexp,
                        FMPRB_RAD_PREC, FMPR_RND_UP);
                    fmpr_add(fmprb_radref(z + xp + yp + k),
                        fmprb_radref(z + xp + yp + k), t,
                        FMPRB_RAD_PREC, FMPR_RND_UP);
                }
            }
        }
    }

    fmpz_clear(zexp);
    fmpr_clear(t);
}

static __inline__ void
_fmprb_poly_addmullow_block(fmprb_ptr z, fmpz * zz,
    const fmpz * xz, const fmpz * xexps, const long * xblocks, long xlen,
    const fmpz * yz, const fmpz * yexps, const long * yblocks, long ylen,
    long n, long prec, int squaring)
{
    long i, j, k, xp, yp, xl, yl, bn;
    fmpz_t zexp;
    fmpr_t t;

    fmpz_init(zexp);
    fmpr_init(t);

    if (squaring)
    {
        for (i = 0; (xp = xblocks[i]) != xlen; i++)
        {
            if (2 * xp >= n)
                continue;

            xl = xblocks[i + 1] - xp;
            bn = FLINT_MIN(2 * xl - 1, n - 2 * xp);
            xl = FLINT_MIN(xl, bn);

            _fmpz_poly_sqrlow(zz, xz + xp, xl, bn);
            fmpz_add_inline(zexp, xexps + i, xexps + i);

            for (k = 0; k < bn; k++)
            {
                fmpr_set_fmpz_2exp(t, zz + k, zexp);
                fmprb_add_fmpr(z + 2 * xp + k, z + 2 * xp + k, t, prec);
            }
        }
    }

    for (i = 0; (xp = xblocks[i]) != xlen; i++)
    {
        for (j = squaring ? i + 1 : 0; (yp = yblocks[j]) != ylen; j++)
        {
            if (xp + yp >= n)
                continue;

            xl = xblocks[i + 1] - xp;
            yl = yblocks[j + 1] - yp;
            bn = FLINT_MIN(xl + yl - 1, n - xp - yp);
            xl = FLINT_MIN(xl, bn);
            yl = FLINT_MIN(yl, bn);

            if (xl >= yl)
                _fmpz_poly_mullow(zz, xz + xp, xl, yz + yp, yl, bn);
            else
                _fmpz_poly_mullow(zz, yz + yp, yl, xz + xp, xl, bn);

           fmpz_add2_fmpz_si_inline(zexp, xexps + i, yexps + j, squaring);

            for (k = 0; k < bn; k++)
            {
                fmpr_set_fmpz_2exp(t, zz + k, zexp);
                fmprb_add_fmpr(z + xp + yp + k, z + xp + yp + k, t, prec);
            }
        }
    }

    fmpz_clear(zexp);
    fmpr_clear(t);
}

void
_fmprb_poly_mullow_block2(fmprb_ptr z, fmprb_srcptr x, long xlen,
                                fmprb_srcptr y, long ylen, long n, long prec)
{
    long xmlen, xrlen, ymlen, yrlen, i;
    fmpz *xz, *yz, *zz;
    fmpz *xe, *ye;
    long *xblocks, *yblocks;
    int squaring;
    fmpz_t scale, t;

    xlen = FLINT_MIN(xlen, n);
    ylen = FLINT_MIN(ylen, n);

    squaring = (x == y) && (xlen == ylen);

    /* Strip trailing zeros */
    xmlen = xrlen = xlen;
    while (xmlen > 0 && fmpr_is_zero(fmprb_midref(x + xmlen - 1))) xmlen--;
    while (xrlen > 0 && fmpr_is_zero(fmprb_radref(x + xrlen - 1))) xrlen--;

    if (squaring)
    {
        ymlen = xmlen;
        yrlen = xrlen;
    }
    else
    {
        ymlen = yrlen = ylen;
        while (ymlen > 0 && fmpr_is_zero(fmprb_midref(y + ymlen - 1))) ymlen--;
        while (yrlen > 0 && fmpr_is_zero(fmprb_radref(y + yrlen - 1))) yrlen--;
    }

    /* We don't know how to deal with infinities or NaNs */
    if (!_fmprb_vec_is_finite(x, xlen) ||
        (!squaring && !_fmprb_vec_is_finite(y, ylen)))
    {
        _fmprb_poly_mullow_classical(z, x, xlen, y, ylen, n, prec);
        return;
    }

    xlen = FLINT_MAX(xmlen, xrlen);
    ylen = FLINT_MAX(ymlen, yrlen);

    /* Start with the zero polynomial */
    _fmprb_vec_zero(z, n);

    /* Nothing to do */
    if (xlen == 0 || ylen == 0)
        return;

    n = FLINT_MIN(n, xlen + ylen - 1);

    fmpz_init(scale);
    fmpz_init(t);
    xz = _fmpz_vec_init(xlen);
    yz = _fmpz_vec_init(ylen);
    zz = _fmpz_vec_init(n);
    xe = _fmpz_vec_init(xlen);
    ye = _fmpz_vec_init(ylen);
    xblocks = flint_malloc(sizeof(long) * (xlen + 1));
    yblocks = flint_malloc(sizeof(long) * (ylen + 1));

    _fmprb_poly_get_scale(scale, x, xlen, y, ylen);

    /* Error propagation */
    /* (xm + xr)*(ym + yr) = (xm*ym) + (xr*ym + xm*yr + xr*yr)
                           = (xm*ym) + (xm*yr + xr*(ym + yr))  */
    if (xrlen != 0 || yrlen != 0)
    {
        fmpr_ptr tmp;
        double *xdbl, *ydbl;
        int can_use_doubles = 1;

        tmp = _fmpr_vec_init(FLINT_MAX(xlen, ylen));
        xdbl = flint_malloc(sizeof(double) * xlen);
        ydbl = flint_malloc(sizeof(double) * ylen);

        /* (xm + xr)^2 = (xm*ym) + (xr^2 + 2 xm xr)
                       = (xm*ym) + xr*(2 xm + xr)    */
        if (squaring)
        {
            can_use_doubles = _fmpr_vec_get_fmpz_2exp_blocks(xz, xdbl, xe,
                xblocks, scale, fmprb_radref(x), xrlen, 2,
                FMPRB_RAD_PREC, can_use_doubles);

            for (i = 0; i < xlen; i++)
            {
                fmpr_abs_round(tmp + i, fmprb_midref(x + i),
                    FMPRB_RAD_PREC, FMPR_RND_UP);
                fmpr_mul_2exp_si(tmp + i, tmp + i, 1);
                fmpr_add(tmp + i, tmp + i, fmprb_radref(x + i),
                    FMPRB_RAD_PREC, FMPR_RND_UP);
            }

            can_use_doubles = _fmpr_vec_get_fmpz_2exp_blocks(yz, ydbl, ye,
                yblocks, scale, tmp, xlen, 1,
                FMPRB_RAD_PREC, can_use_doubles);

            _fmprb_poly_addmullow_rad(z, zz,
                xz, xdbl, xe, xblocks, xrlen,
                yz, ydbl, ye, yblocks, xlen, n, can_use_doubles);
        }
        else if (yrlen == 0)
        {
            /* xr * |ym| */
            can_use_doubles = _fmpr_vec_get_fmpz_2exp_blocks(xz, xdbl, xe,
                xblocks, scale, fmprb_radref(x), xrlen, 2,
                FMPRB_RAD_PREC, can_use_doubles);

            for (i = 0; i < ymlen; i++)
                fmpr_abs_round(tmp + i, fmprb_midref(y + i),
                    FMPRB_RAD_PREC, FMPR_RND_UP);

            can_use_doubles = _fmpr_vec_get_fmpz_2exp_blocks(yz, ydbl, ye,
                yblocks, scale, tmp, ymlen, 1,
                FMPRB_RAD_PREC, can_use_doubles);

            _fmprb_poly_addmullow_rad(z, zz,
                xz, xdbl, xe, xblocks, xrlen,
                yz, ydbl, ye, yblocks, ymlen, n, can_use_doubles);
        }
        else
        {
            /* |xm| * yr */
            for (i = 0; i < xmlen; i++)
                fmpr_abs_round(tmp + i, fmprb_midref(x + i),
                    FMPRB_RAD_PREC, FMPR_RND_UP);

            can_use_doubles = _fmpr_vec_get_fmpz_2exp_blocks(xz, xdbl, xe,
                xblocks, scale, tmp, xmlen, 1,
                FMPRB_RAD_PREC, can_use_doubles);

            can_use_doubles = _fmpr_vec_get_fmpz_2exp_blocks(yz, ydbl, ye,
                yblocks, scale, fmprb_radref(y), yrlen, 2,
                FMPRB_RAD_PREC, can_use_doubles);

            _fmprb_poly_addmullow_rad(z, zz,
                xz, xdbl, xe, xblocks, xmlen,
                yz, ydbl, ye, yblocks, yrlen, n, can_use_doubles);

            /* xr*(|ym| + yr) */
            if (xrlen != 0)
            {
                can_use_doubles = 1;
                can_use_doubles = _fmpr_vec_get_fmpz_2exp_blocks(xz, xdbl, xe,
                    xblocks, scale, fmprb_radref(x), xrlen, 2,
                    FMPRB_RAD_PREC, can_use_doubles);

                for (i = 0; i < ylen; i++)
                    fmpr_add_abs_ubound(tmp + i, fmprb_midref(y + i),
                        fmprb_radref(y + i), FMPRB_RAD_PREC);

                can_use_doubles = _fmpr_vec_get_fmpz_2exp_blocks(yz, ydbl, ye,
                    yblocks, scale, tmp, ylen, 1,
                    FMPRB_RAD_PREC, can_use_doubles);

                _fmprb_poly_addmullow_rad(z, zz,
                    xz, xdbl, xe, xblocks, xrlen,
                    yz, ydbl, ye, yblocks, ylen, n, can_use_doubles);
            }
        }

        _fmpr_vec_clear(tmp, FLINT_MAX(xlen, ylen));
        flint_free(xdbl);
        flint_free(ydbl);
    }

    /* multiply midpoints */
    if (xmlen != 0 && ymlen != 0)
    {
        _fmpr_vec_get_fmpz_2exp_blocks(xz, NULL, xe, xblocks,
            scale, fmprb_midref(x), xmlen, 2, prec, 0);

        if (squaring)
        {
            _fmprb_poly_addmullow_block(z, zz,
                xz, xe, xblocks, xmlen, xz, xe, xblocks, xmlen, n, prec, 1);
        }
        else
        {
            _fmpr_vec_get_fmpz_2exp_blocks(yz, NULL, ye, yblocks,
                scale, fmprb_midref(y), ymlen, 2, prec, 0);

            _fmprb_poly_addmullow_block(z, zz,
                xz, xe, xblocks, xmlen,
                yz, ye, yblocks, ymlen, n, prec, 0);
        }
    }

    /* Unscale. */
    if (!fmpz_is_zero(scale))
    {
        fmpz_zero(t);
        for (i = 0; i < n; i++)
        {
            fmprb_mul_2exp_fmpz(z + i, z + i, t);
            fmpz_add(t, t, scale);
        }
    }

    _fmpz_vec_clear(xz, xlen);
    _fmpz_vec_clear(yz, ylen);
    _fmpz_vec_clear(zz, n);
    _fmpz_vec_clear(xe, xlen);
    _fmpz_vec_clear(ye, ylen);
    flint_free(xblocks);
    flint_free(yblocks);
    fmpz_clear(scale);
    fmpz_clear(t);
}

void
fmprb_poly_mullow_block2(fmprb_poly_t res, const fmprb_poly_t poly1,
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
        _fmprb_poly_mullow_block2(tmp->coeffs, poly1->coeffs, xlen,
            poly2->coeffs, ylen, zlen, prec);
        fmprb_poly_swap(res, tmp);
        fmprb_poly_clear(tmp);
    }
    else
    {
        fmprb_poly_fit_length(res, zlen);
        _fmprb_poly_mullow_block2(res->coeffs, poly1->coeffs, xlen,
            poly2->coeffs, ylen, zlen, prec);
    }

    _fmprb_poly_set_length(res, zlen);
    _fmprb_poly_normalise(res);
}

