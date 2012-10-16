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

#include "fmpcb_poly.h"

/* we don't need any error bounding, so we define a few helper
   functions that ignore the radii */

static __inline__ void
fmpcb_sub_mid(fmpcb_t z, const fmpcb_t x, const fmpcb_t y, long prec)
{
    fmpr_sub(fmprb_midref(fmpcb_realref(z)),
        fmprb_midref(fmpcb_realref(x)),
        fmprb_midref(fmpcb_realref(y)), prec, FMPR_RND_DOWN);
    fmpr_sub(fmprb_midref(fmpcb_imagref(z)),
        fmprb_midref(fmpcb_imagref(x)),
        fmprb_midref(fmpcb_imagref(y)), prec, FMPR_RND_DOWN);
}

static __inline__ void
fmpcb_add_mid(fmpcb_t z, const fmpcb_t x, const fmpcb_t y, long prec)
{
    fmpr_add(fmprb_midref(fmpcb_realref(z)),
        fmprb_midref(fmpcb_realref(x)),
        fmprb_midref(fmpcb_realref(y)), prec, FMPR_RND_DOWN);
    fmpr_add(fmprb_midref(fmpcb_imagref(z)),
        fmprb_midref(fmpcb_imagref(x)),
        fmprb_midref(fmpcb_imagref(y)), prec, FMPR_RND_DOWN);
}

static __inline__ void
fmpcb_mul_mid(fmpcb_t z, const fmpcb_t x, const fmpcb_t y, long prec)
{
#define a fmprb_midref(fmpcb_realref(x))
#define b fmprb_midref(fmpcb_imagref(x))
#define c fmprb_midref(fmpcb_realref(y))
#define d fmprb_midref(fmpcb_imagref(y))
#define e fmprb_midref(fmpcb_realref(z))
#define f fmprb_midref(fmpcb_imagref(z))

    fmpr_t t, u, v;

    fmpr_init(t);
    fmpr_init(u);
    fmpr_init(v);

    fmpr_add(t, a, b, prec, FMPR_RND_DOWN);
    fmpr_add(u, c, d, prec, FMPR_RND_DOWN);
    fmpr_mul(v, t, u, prec, FMPR_RND_DOWN);

    fmpr_mul(t, a, c, prec, FMPR_RND_DOWN);
    fmpr_mul(u, b, d, prec, FMPR_RND_DOWN);

    fmpr_sub(e, t, u, prec, FMPR_RND_DOWN);
    fmpr_sub(f, v, t, prec, FMPR_RND_DOWN);
    fmpr_sub(f, f, u, prec, FMPR_RND_DOWN);

    fmpr_clear(t);
    fmpr_clear(u);
    fmpr_clear(v);

#undef a
#undef b
#undef c
#undef d
#undef e
#undef f
}

static __inline__ void
fmpcb_inv_mid(fmpcb_t z, const fmpcb_t x, long prec)
{
    fmpr_t t;
    fmpr_init(t);

#define a fmprb_midref(fmpcb_realref(x))
#define b fmprb_midref(fmpcb_imagref(x))
#define e fmprb_midref(fmpcb_realref(z))
#define f fmprb_midref(fmpcb_imagref(z))

    fmpr_mul(t, a, a, prec, FMPR_RND_DOWN);
    fmpr_addmul(t, b, b, prec, FMPR_RND_DOWN);

    fmpr_div(e, a, t, prec, FMPR_RND_DOWN);
    fmpr_div(f, b, t, prec, FMPR_RND_DOWN);

    fmpr_neg(f, f);

#undef a
#undef b
#undef e
#undef f

    fmpr_clear(t);
}

void
_fmpcb_poly_evaluate_mid(fmpcb_t res, const fmpcb_struct * f, long len,
                           const fmpcb_t a, long prec)
{
    long i = len - 1;
    fmpcb_t t;

    fmpcb_init(t);
    fmpcb_set(res, f + i);

    for (i = len - 2; i >= 0; i--)
    {
        fmpcb_mul_mid(t, res, a, prec);
        fmpcb_add_mid(res, f + i, t, prec);
    }

    fmpcb_clear(t);
}

void
_fmpcb_poly_refine_roots_durand_kerner(fmpcb_struct * roots,
        const fmpcb_struct * poly, long len, long prec)
{
    long i, j;

    fmpcb_t x, y, t;

    fmpcb_init(x);
    fmpcb_init(y);
    fmpcb_init(t);

    for (i = 0; i < len - 1; i++)
    {
        _fmpcb_poly_evaluate_mid(x, poly, len, roots + i, prec);

        fmpcb_set(y, poly + len - 1);

        for (j = 0; j < len - 1; j++)
        {
            if (i != j)
            {
                fmpcb_sub_mid(t, roots + i, roots + j, prec);
                fmpcb_mul_mid(y, y, t, prec);
            }
        }

        fmpr_zero(fmprb_radref(fmpcb_realref(y)));
        fmpr_zero(fmprb_radref(fmpcb_imagref(y)));

        fmpcb_inv_mid(t, y, prec);
        fmpcb_mul_mid(t, t, x, prec);

        fmpcb_sub_mid(roots + i, roots + i, t, prec);

        fmpr_set_round(fmprb_radref(fmpcb_realref(roots + i)),
            fmprb_midref(fmpcb_realref(t)), FMPRB_RAD_PREC, FMPR_RND_UP);
        fmpr_set_round(fmprb_radref(fmpcb_imagref(roots + i)),
            fmprb_midref(fmpcb_imagref(t)), FMPRB_RAD_PREC, FMPR_RND_UP);
    }

    fmpcb_clear(x);
    fmpcb_clear(y);
    fmpcb_clear(t);
}

