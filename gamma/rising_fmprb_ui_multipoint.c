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

    Copyright (C) 2012, 2013 Fredrik Johansson

******************************************************************************/

#include "gamma.h"
#include "fmprb_poly.h"

void
gamma_rising_fmprb_ui_multipoint(fmprb_t f, const fmprb_t c, ulong n, long prec)
{
    long i, m, wp;
    fmprb_struct *t, *u, *v;
    fmprb_t r, w;
    ulong s;

    if (n <= 1)
    {
        if (n == 0)
            fmprb_one(f);
        else
            fmprb_set_round(f, c, prec);
        return;
    }

    /* TODO: this is useless if the input isn't as precise as the working
       precision to begin with. We rather want to convert to an exact input
       within the evaluation, and add the propagated error afterwards. */
    wp = FMPR_PREC_ADD(prec, n);

    fmprb_init(r);
    fmprb_init(w);

    m = n_sqrt(n);

    t = _fmprb_vec_init(m + 1);
    u = _fmprb_vec_init(m + 1);
    v = _fmprb_vec_init(m + 1);

    /* the polynomial is x(x+1)(x+2)... */
    for (i = 0; i < m; i++)
        fmprb_set_si(t + i, -i);

    _fmprb_poly_product_roots(u, t, m, wp);

    /* the evaluation points are c, c+m, ... */
    for (i = 0; i < m; i++)
        fmprb_add_ui(t + i, c, i * m, wp);

    _fmprb_poly_evaluate_vec_fast(v, u, m + 1, t, m, wp);

    fmprb_one(r);
    for (i = 0; i < m; i++)
        fmprb_mul(r, r, v + i, wp);

    /* remaining part of product (todo: this should use binary
       splitting as well, to minimize numerical error) */
    for (s = m * m; s < n; s++)
    {
        fmprb_add_ui(w, c, s, wp);
        fmprb_mul(r, r, w, wp);
    }

    fmprb_set(f, r);

    _fmprb_vec_clear(t, m + 1);
    _fmprb_vec_clear(u, m + 1);
    _fmprb_vec_clear(v, m + 1);

    fmprb_clear(r);
    fmprb_clear(w);
}

