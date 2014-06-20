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

#include "arf.h"

int _arf_are_close(const arf_t x, const arf_t y, long prec)
{
    fmpz_t xb, yb;
    fmpz_t delta;
    int result;

    fmpz_init(xb);
    fmpz_init(yb);
    fmpz_init(delta);

    arf_bot(xb, x);
    arf_bot(yb, y);

    if (fmpz_cmp(ARF_EXPREF(x), ARF_EXPREF(y)) >= 0)
        fmpz_sub(delta, xb, ARF_EXPREF(y));
    else
        fmpz_sub(delta, yb, ARF_EXPREF(x));

    fmpz_sub_ui(delta, delta, 64);
    result = (fmpz_cmp_ui(delta, prec) < 0);

    fmpz_clear(xb);
    fmpz_clear(yb);
    fmpz_clear(delta);

    return result;
}

int
_arf_add_eps(arf_t s, const arf_t x, int sgn, long prec, arf_rnd_t rnd)
{
    arf_t t;
    long bits;

    bits = arf_bits(x);

    if (bits == 0)
    {
        printf("_arf_add_eps\n");
        abort();
    }

    bits = FLINT_MAX(bits, prec) + 10;

    arf_init(t);
    arf_set_si(t, sgn);
    arf_mul_2exp_fmpz(t, t, ARF_EXPREF(x));
    arf_mul_2exp_si(t, t, -bits);
    arf_add(s, x, t, prec, rnd);
    arf_clear(t);

    return 1;
}

int
arf_sum(arf_t s, arf_srcptr terms, long len, long prec, arf_rnd_t rnd)
{
    arf_ptr blocks;
    long i, j, used;
    int have_merged, res;

    /* first check if the result is inf or nan */
    {
        int have_pos_inf = 0;
        int have_neg_inf = 0;

        for (i = 0; i < len; i++)
        {
            if (arf_is_pos_inf(terms + i))
            {
                if (have_neg_inf)
                {
                    arf_nan(s);
                    return 0;
                }
                have_pos_inf = 1;
            }
            else if (arf_is_neg_inf(terms + i))
            {
                if (have_pos_inf)
                {
                    arf_nan(s);
                    return 0;
                }
                have_neg_inf = 1;
            }
            else if (arf_is_nan(terms + i))
            {
                arf_nan(s);
                return 0;
            }
        }

        if (have_pos_inf)
        {
            arf_pos_inf(s);
            return 0;
        }

        if (have_neg_inf)
        {
            arf_neg_inf(s);
            return 0;
        }
    }

    blocks = flint_malloc(sizeof(arf_struct) * len);
    for (i = 0; i < len; i++)
        arf_init(blocks + i);

    /* put all terms into blocks */
    used = 0;
    for (i = 0; i < len; i++)
    {
        if (!arf_is_zero(terms + i))
        {
            arf_set(blocks + used, terms + i);
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
                if (_arf_are_close(blocks + i, blocks + j, prec))
                {
                    arf_add(blocks + i, blocks + i, blocks + j,
                        ARF_PREC_EXACT, ARF_RND_DOWN);

                    /* remove the merged block */
                    arf_swap(blocks + j, blocks + used - 1);
                    used--;

                    /* remove the updated block if the sum is zero */
                    if (arf_is_zero(blocks + i))
                    {
                        arf_swap(blocks + i, blocks + used - 1);
                        used--;
                    }

                    have_merged = 1;
                }
            }
        }
    }

    if (used == 0)
    {
        arf_zero(s);
        res = 0;
    }
    else if (used == 1)
    {
        res = arf_set_round(s, blocks + 0, prec, rnd);
    }
    else
    {
        /* find the two largest blocks */
        for (i = 1; i < used; i++)
            if (arf_cmpabs(blocks + 0, blocks + i) < 0)
                arf_swap(blocks + 0, blocks + i);

        for (i = 2; i < used; i++)
            if (arf_cmpabs(blocks + 1, blocks + i) < 0)
                arf_swap(blocks + 1, blocks + i);

        res = _arf_add_eps(s, blocks + 0, arf_sgn(blocks + 1), prec, rnd);
    }

    for (i = 0; i < len; i++)
        arf_clear(blocks + i);
    flint_free(blocks);

    return res;
}

