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

#include "acb_modular.h"

int
acb_modular_is_in_fundamental_domain(const acb_t z, const arf_t tol, long prec)
{
    arb_t t;
    arb_init(t);

    /* require re(w) + 1/2 >= 0 */
    arb_set_ui(t, 1);
    arb_mul_2exp_si(t, t, -1);
    arb_add(t, t, acb_realref(z), prec);
    arb_add_arf(t, t, tol, prec);
    if (!arb_is_nonnegative(t))
    {
        arb_clear(t);
        return 0;
    }

    /* require re(w) - 1/2 <= 0 */
    arb_set_ui(t, 1);
    arb_mul_2exp_si(t, t, -1);
    arb_sub(t, acb_realref(z), t, prec);
    arb_sub_arf(t, t, tol, prec);
    if (!arb_is_nonpositive(t))
    {
        arb_clear(t);
        return 0;
    }

    /* require |w| >= 1 - tol, i.e. |w| - 1 + tol >= 0 */
    acb_abs(t, z, prec);
    arb_sub_ui(t, t, 1, prec);
    arb_add_arf(t, t, tol, prec);
    if (!arb_is_nonnegative(t))
    {
        arb_clear(t);
        return 0;
    }

    arb_clear(t);
    return 1;
}

