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

#include "fmpr.h"

int are_close(const fmpr_t x, const fmpr_t y, long prec)
{
    fmpz_t xb, yb;
    fmpz_t delta;
    int result;

    fmpz_init(xb);
    fmpz_init(yb);
    fmpz_init(delta);

    fmpz_add_ui(xb, fmpr_expref(x), fmpz_bits(fmpr_manref(x)));
    fmpz_add_ui(yb, fmpr_expref(y), fmpz_bits(fmpr_manref(y)));

    if (fmpz_cmp(xb, yb) >= 0)
        fmpz_sub(delta, fmpr_expref(x), yb);
    else
        fmpz_sub(delta, fmpr_expref(y), xb);

    fmpz_sub_ui(delta, delta, 64);
    result = (fmpz_cmp_ui(delta, prec) < 0);

    fmpz_clear(xb);
    fmpz_clear(yb);
    fmpz_clear(delta);

    return result;
}

long
fmpr_sum(fmpr_t s, const fmpr_struct * terms, long len, long prec, fmpr_rnd_t rnd)
{
    fmpr_struct * blocks;
    long i, j, used, res;
    int have_merged;

    /* first check if the result is inf or nan */
    {
        int have_pos_inf = 0;
        int have_neg_inf = 0;

        for (i = 0; i < len; i++)
        {
            if (fmpr_is_pos_inf(terms + i))
            {
                if (have_neg_inf)
                {
                    fmpr_nan(s);
                    return FMPR_RESULT_EXACT;
                }
                have_pos_inf = 1;
            }
            else if (fmpr_is_neg_inf(terms + i))
            {
                if (have_pos_inf)
                {
                    fmpr_nan(s);
                    return FMPR_RESULT_EXACT;
                }
                have_neg_inf = 1;
            }
            else if (fmpr_is_nan(terms + i))
            {
                fmpr_nan(s);
                return FMPR_RESULT_EXACT;
            }
        }

        if (have_pos_inf)
        {
            fmpr_pos_inf(s);
            return FMPR_RESULT_EXACT;
        }

        if (have_neg_inf)
        {
            fmpr_neg_inf(s);
            return FMPR_RESULT_EXACT;
        }
    }

    blocks = flint_malloc(sizeof(fmpr_struct) * len);
    for (i = 0; i < len; i++)
        fmpr_init(blocks + i);

    /* put all terms into blocks */
    used = 0;
    for (i = 0; i < len; i++)
    {
        if (!fmpr_is_zero(terms + i))
        {
            fmpr_set(blocks + used, terms + i);
            used++;
        }
    }

    /* merge blocks until all are well separated */
    have_merged = 1;
    while (used >= 2 && have_merged)
    {
        have_merged = 0;

        for (i = 0; i < used && !have_merged; i++)
        {
            for (j = i + 1; j < used && !have_merged; j++)
            {
                if (are_close(blocks + i, blocks + j, prec))
                {
                    fmpr_add(blocks + i, blocks + i, blocks + j, FMPR_PREC_EXACT, FMPR_RND_DOWN);

                    /* remove the merged block */
                    fmpr_swap(blocks + j, blocks + used - 1);
                    used--;

                    /* remove the updated block if the sum is zero */
                    if (fmpr_is_zero(blocks + i))
                    {
                        fmpr_swap(blocks + i, blocks + used - 1);
                        used--;
                    }

                    have_merged = 1;
                }
            }
        }
    }

    if (used == 0)
    {
        fmpr_zero(s);
        res = FMPR_RESULT_EXACT;
    }
    else if (used == 1)
    {
        res = fmpr_set_round(s, blocks + 0, prec, rnd);
    }
    else
    {
        /* find the two largest blocks */
        for (i = 1; i < used; i++)
            if (fmpr_cmpabs(blocks + 0, blocks + i) < 0)
                fmpr_swap(blocks + 0, blocks + i);

        for (i = 2; i < used; i++)
            if (fmpr_cmpabs(blocks + 1, blocks + i) < 0)
                fmpr_swap(blocks + 1, blocks + i);

        res = _fmpr_add_eps(s, blocks + 0, fmpr_sgn(blocks + 1), prec, rnd);
    }

    for (i = 0; i < len; i++)
        fmpr_clear(blocks + i);
    flint_free(blocks);

    return res;
}

