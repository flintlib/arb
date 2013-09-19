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

long
_fmpcb_get_mid_mag(const fmpcb_t z)
{
    long rm, im;

    rm = fmpr_abs_bound_lt_2exp_si(fmprb_midref(fmpcb_realref(z)));
    im = fmpr_abs_bound_lt_2exp_si(fmprb_midref(fmpcb_imagref(z)));

    return FLINT_MAX(rm, im);
}

long
_fmpcb_get_rad_mag(const fmpcb_t z)
{
    long rm, im;

    rm = fmpr_abs_bound_lt_2exp_si(fmprb_radref(fmpcb_realref(z)));
    im = fmpr_abs_bound_lt_2exp_si(fmprb_radref(fmpcb_imagref(z)));

    return FLINT_MAX(rm, im);
}

void
_fmpcb_poly_roots_initial_values(fmpcb_ptr roots, long deg, long prec)
{
    long i;

    fmpq_t q;
    fmpq_init(q);

    fmpq_set_si(q, 4, 10);
    fmprb_set_fmpq(fmpcb_realref(roots + 0), q, prec);
    fmpr_zero(fmprb_radref(fmpcb_realref(roots + 0)));

    fmpq_set_si(q, 9, 10);
    fmprb_set_fmpq(fmpcb_imagref(roots + 0), q, prec);
    fmpr_zero(fmprb_radref(fmpcb_imagref(roots + 0)));
    fmpq_clear(q);

    for (i = 1; i < deg; i++)
    {
        fmpcb_mul(roots + i, roots + i - 1, roots + 0, prec);
        fmpr_zero(fmprb_radref(fmpcb_realref(roots + i)));
        fmpr_zero(fmprb_radref(fmpcb_imagref(roots + i)));
    }
}

long
_fmpcb_poly_find_roots(fmpcb_ptr roots,
    fmpcb_srcptr poly,
    fmpcb_srcptr initial, long len, long maxiter, long prec)
{
    long iter, i, deg;
    long rootmag, max_rootmag, correction, max_correction;

    deg = len - 1;

    if (deg == 0)
    {
        return 0;
    }
    else if (fmpcb_contains_zero(poly + len - 1))
    {
        /* if the leading coefficient contains zero, roots can be anywhere */
        for (i = 0; i < deg; i++)
        {
            fmprb_zero(fmpcb_realref(roots + i));
            fmpr_pos_inf(fmprb_radref(fmpcb_realref(roots + i)));
            fmprb_zero(fmpcb_imagref(roots + i));
            fmpr_pos_inf(fmprb_radref(fmpcb_imagref(roots + i)));
        }
        return 0;
    }
    else if (deg == 1)
    {
        fmpcb_inv(roots + 0, poly + 1, prec);
        fmpcb_mul(roots + 0, roots + 0, poly + 0, prec);
        fmpcb_neg(roots + 0, roots + 0);
        return 1;
    }

    if (initial == NULL)
        _fmpcb_poly_roots_initial_values(roots, deg, prec);
    else
        _fmpcb_vec_set(roots, initial, deg);

    if (maxiter == 0)
        maxiter = 2 * deg + n_sqrt(prec);

    for (iter = 0; iter < maxiter; iter++)
    {
        max_rootmag = -FMPR_PREC_EXACT;
        for (i = 0; i < deg; i++)
        {
            rootmag = _fmpcb_get_mid_mag(roots + i);
            max_rootmag = FLINT_MAX(rootmag, max_rootmag);
        }

        _fmpcb_poly_refine_roots_durand_kerner(roots, poly, len, prec);

        max_correction = -FMPR_PREC_EXACT;
        for (i = 0; i < deg; i++)
        {
            correction = _fmpcb_get_rad_mag(roots + i);
            max_correction = FLINT_MAX(correction, max_correction);
        }

        /* estimate the correction relative to the whole set of roots */
        max_correction -= max_rootmag;

        /* printf("ITER %ld MAX CORRECTION: %ld\n", iter, max_correction); */

        if (max_correction < -prec / 2)
            maxiter = FLINT_MIN(maxiter, iter + 2);
        else if (max_correction < -prec / 3)
            maxiter = FLINT_MIN(maxiter, iter + 3);
        else if (max_correction < -prec / 4)
            maxiter = FLINT_MIN(maxiter, iter + 4);
    }

    return _fmpcb_poly_validate_roots(roots, poly, len, prec);
}


long
fmpcb_poly_find_roots(fmpcb_ptr roots,
    const fmpcb_poly_t poly, fmpcb_srcptr initial,
    long maxiter, long prec)
{
    long len = poly->length;

    if (len == 0)
    {
        printf("find_roots: expected a nonzero polynomial");
        abort();
    }

    return _fmpcb_poly_find_roots(roots, poly->coeffs, initial,
                len, maxiter, prec);
}

