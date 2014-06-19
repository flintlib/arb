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

#include "arb_calc.h"

#define BLOCK_NO_ZERO 0
#define BLOCK_ISOLATED_ZERO 1
#define BLOCK_UNKNOWN 2

/* 0 means that it *could* be zero; otherwise +/- 1 */
static __inline__ int
_arb_sign(const arb_t t)
{
    if (arb_is_positive(t))
        return 1;
    else if (arb_is_negative(t))
        return -1;
    else
        return 0;
}

static int
check_block(arb_calc_func_t func, void * param, const arf_interval_t block,
    int asign, int bsign, long prec)
{
    arb_struct t[2];
    arb_t x;
    int result;

    arb_init(t + 0);
    arb_init(t + 1);
    arb_init(x);

    arf_interval_get_arb(x, block, prec);
    func(t, x, param, 1, prec);

    result = BLOCK_UNKNOWN;

    if (arb_is_positive(t) || arb_is_negative(t))
    {
        result = BLOCK_NO_ZERO;
    }
    else
    {
        if ((asign < 0 && bsign > 0) || (asign > 0 && bsign < 0))
        {
            func(t, x, param, 2, prec);

            if (arb_is_finite(t + 1) && !arb_contains_zero(t + 1))
            {
                result = BLOCK_ISOLATED_ZERO;
            }
        }
    }

    arb_clear(t + 0);
    arb_clear(t + 1);
    arb_clear(x);

    return result;
}

#define ADD_BLOCK       \
    if (*length >= *alloc)   \
    {   \
        long new_alloc;   \
        new_alloc = (*alloc == 0) ? 1 : 2 * (*alloc); \
        *blocks = flint_realloc(*blocks, sizeof(arf_interval_struct) * new_alloc);   \
        *flags = flint_realloc(*flags, sizeof(int) * new_alloc);   \
        *alloc = new_alloc;   \
    }   \
    arf_interval_init((*blocks) + *length);   \
    arf_interval_set((*blocks) + *length, block);   \
    (*flags)[*length] = status; \
    (*length)++; \


static void
isolate_roots_recursive(arf_interval_ptr * blocks, int ** flags,
    long * length, long * alloc,
    arb_calc_func_t func, void * param,
    const arf_interval_t block, int asign, int bsign,
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
                    if (arb_calc_verbose)
                    {
                        printf("found isolated root in: ");
                        arf_interval_printd(block, 15);
                        printf("\n");
                    }

                    *found_count -= 1;
                }

                ADD_BLOCK
            }
            else
            {
                arf_interval_t L, R;
                int msign;

                arf_interval_init(L);
                arf_interval_init(R);

                msign = arb_calc_partition(L, R, func, param, block, prec);

                if (msign == 0 && arb_calc_verbose)
                {
                    printf("possible zero at midpoint: ");
                    arf_interval_printd(block, 15);
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

                arf_interval_clear(L);
                arf_interval_clear(R);
            }
        }
    }
}

long
arb_calc_isolate_roots(arf_interval_ptr * blocks, int ** flags,
    arb_calc_func_t func, void * param,
    const arf_interval_t block, long maxdepth, long maxeval, long maxfound,
    long prec)
{
    int asign, bsign;
    long length, alloc;
    arb_t m, v;

    *blocks = NULL;
    *flags = NULL;
    length = 0;
    alloc = 0;

    arb_init(m);
    arb_init(v);

    arb_set_arf(m, &block->a);
    func(v, m, param, 1, prec);
    asign = _arb_sign(v);

    arb_set_arf(m, &block->b);
    func(v, m, param, 1, prec);
    bsign = _arb_sign(v);

    arb_clear(m);
    arb_clear(v);

    isolate_roots_recursive(blocks, flags, &length, &alloc,
        func, param, block, asign, bsign,
        maxdepth, &maxeval, &maxfound, prec);

    *blocks = flint_realloc(*blocks, length * sizeof(arf_interval_struct));
    *flags = flint_realloc(*flags, length * sizeof(int));

    return length;
}

