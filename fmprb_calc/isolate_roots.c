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

#include "fmprb_calc.h"

#define BLOCK_NO_ZERO 0
#define BLOCK_ISOLATED_ZERO 1
#define BLOCK_UNKNOWN 2

/* 0 means that it *could* be zero; otherwise +/- 1 */
static __inline__ int
_fmprb_sign(const fmprb_t t)
{
    if (fmprb_is_positive(t))
        return 1;
    else if (fmprb_is_negative(t))
        return -1;
    else
        return 0;
}

static int
check_block(fmprb_calc_func_t func, void * param, const fmprb_t block,
    int asign, int bsign, long prec)
{
    fmprb_struct t[2];
    int result;

    fmprb_init(t + 0);
    fmprb_init(t + 1);

    func(t, block, param, 1, prec);

    result = BLOCK_UNKNOWN;

    if (fmprb_is_positive(t) || fmprb_is_negative(t))
    {
        result = BLOCK_NO_ZERO;
    }
    else
    {
        if ((asign < 0 && bsign > 0) || (asign > 0 && bsign < 0))
        {
            func(t, block, param, 2, prec);

            if (fmprb_is_finite(t + 1) && !fmprb_contains_zero(t + 1))
            {
                result = BLOCK_ISOLATED_ZERO;
            }
        }
    }

    fmprb_clear(t + 0);
    fmprb_clear(t + 1);

    return result;
}

static int
partition(fmprb_t L, fmprb_t R,
    fmprb_calc_func_t func, void * param, const fmprb_t block, long prec)
{
    fmprb_t t, m;
    int msign;

    fmprb_init(t);
    fmprb_init(m);

    fmprb_set_fmpr(m, fmprb_midref(block));
    func(t, m, param, 1, prec);
    msign = _fmprb_sign(t);

    fmpr_mul_2exp_si(fmprb_radref(L), fmprb_radref(block), -1);
    fmpr_set(fmprb_radref(R), fmprb_radref(L));

    /* XXX: deal with huge shifts */
    fmpr_sub(fmprb_midref(L), fmprb_midref(block), fmprb_radref(L), FMPR_PREC_EXACT, FMPR_RND_DOWN);
    fmpr_add(fmprb_midref(R), fmprb_midref(block), fmprb_radref(R), FMPR_PREC_EXACT, FMPR_RND_DOWN);

    fmprb_clear(t);
    fmprb_clear(m);

    return msign;
}

#define ADD_BLOCK       \
    if (*length >= *alloc)   \
    {   \
        long new_alloc;   \
        new_alloc = (*alloc == 0) ? 1 : 2 * (*alloc); \
        *blocks = flint_realloc(*blocks, sizeof(fmprb_struct) * new_alloc);   \
        *flags = flint_realloc(*flags, sizeof(int) * new_alloc);   \
        *alloc = new_alloc;   \
    }   \
    fmprb_init((*blocks) + *length);   \
    fmprb_set((*blocks) + *length, block);   \
    (*flags)[*length] = status; \
    (*length)++; \


static void
isolate_roots_recursive(fmprb_ptr * blocks, int ** flags,
    long * length, long * alloc,
    fmprb_calc_func_t func, void * param,
    const fmprb_t block, int asign, int bsign,
    long depth, long * eval_count, long * found_count,
    long prec)
{
    int status;

    if (*found_count <= 0 || *eval_count <= 0)
    {
        status = BLOCK_UNKNOWN;
        ADD_BLOCK
    }
    else
    {
        *eval_count -= 1;
        status = check_block(func, param, block, asign, bsign, prec);

        if (status != BLOCK_NO_ZERO)
        {
            if (status == BLOCK_ISOLATED_ZERO || depth <= 0)
            {
                if (status == BLOCK_ISOLATED_ZERO)
                {
                    if (fmprb_calc_verbose)
                    {
                        printf("found isolated root in: ");
                        fmprb_printd(block, 15);
                        printf("\n");
                    }

                    *found_count -= 1;
                }

                ADD_BLOCK
            }
            else
            {
                fmprb_t L, R;
                int msign;

                fmprb_init(L);
                fmprb_init(R);

                msign = partition(L, R, func, param, block, prec);

                if (msign == 0 && fmprb_calc_verbose)
                {
                    printf("possible zero at midpoint: ");
                    fmprb_printd(block, 15);
                    printf("\n");
                }

                isolate_roots_recursive(blocks, flags, length, alloc,
                    func, param,
                    L, asign, msign,
                    depth - 1, eval_count, found_count, prec);

                isolate_roots_recursive(blocks, flags, length, alloc,
                    func, param,
                    R, msign, bsign,
                    depth - 1, eval_count, found_count, prec);

                fmprb_clear(L);
                fmprb_clear(R);
            }
        }
    }
}

long
fmprb_calc_isolate_roots(fmprb_ptr * blocks, int ** flags,
    fmprb_calc_func_t func, void * param,
    const fmprb_t block, long maxdepth, long maxeval, long maxfound,
    long prec)
{
    int asign, bsign;
    long length, alloc;
    fmprb_t t, u;

    *blocks = NULL;
    *flags = NULL;
    length = 0;
    alloc = 0;

    fmprb_init(t);
    fmprb_init(u);

    /* XXX: deal with huge shifts */
    fmpr_sub(fmprb_midref(t), fmprb_midref(block), fmprb_radref(block), FMPR_PREC_EXACT, FMPR_RND_DOWN);
    func(u, t, param, 1, prec);
    asign = _fmprb_sign(u);

    fmpr_add(fmprb_midref(t), fmprb_midref(block), fmprb_radref(block), FMPR_PREC_EXACT, FMPR_RND_DOWN);
    func(u, t, param, 1, prec);
    bsign = _fmprb_sign(u);

    fmprb_clear(t);
    fmprb_clear(u);

    isolate_roots_recursive(blocks, flags, &length, &alloc,
        func, param, block, asign, bsign,
        maxdepth, &maxeval, &maxfound, prec);

    *blocks = flint_realloc(*blocks, length * sizeof(fmprb_struct));
    *flags = flint_realloc(*flags, length * sizeof(int));

    return length;
}

