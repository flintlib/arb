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

#include "fmprb.h"

#define ASSERT(cond) if (!(cond)) { printf("FAIL: %d\n", __LINE__); abort(); }

int main()
{
    fmprb_t zero, pos, neg, pos_inf, neg_inf, pos_inf_err, neg_inf_err,
      zero_pm_inf, pos_pm_inf, neg_pm_inf,
      indet_exact, indet_pos_rad, indet_inf_rad;

    printf("special....");
    fflush(stdout);

    fmprb_init(zero);
    fmprb_init(pos);
    fmprb_init(neg);
    fmprb_init(pos_inf);
    fmprb_init(neg_inf);
    fmprb_init(pos_inf_err);
    fmprb_init(neg_inf_err);
    fmprb_init(zero_pm_inf);
    fmprb_init(pos_pm_inf);
    fmprb_init(neg_pm_inf);
    fmprb_init(indet_exact);
    fmprb_init(indet_pos_rad);
    fmprb_init(indet_inf_rad);

    fmprb_set_si(pos, 3); fmprb_div_ui(pos, pos, 5, 53);
    fmprb_neg(neg, pos);
    fmprb_pos_inf(pos_inf);
    fmprb_neg_inf(neg_inf);
    fmprb_pos_inf(pos_inf_err); fmpr_set_ui(fmprb_radref(pos_inf_err), 3);
    fmprb_neg(neg_inf_err, pos_inf_err);
    fmprb_zero_pm_inf(zero_pm_inf);
    fmpr_set_si(fmprb_midref(pos_pm_inf), 3); fmpr_pos_inf(fmprb_radref(pos_pm_inf));
    fmpr_set_si(fmprb_midref(neg_pm_inf), -3); fmpr_pos_inf(fmprb_radref(neg_pm_inf));
    fmpr_nan(fmprb_midref(indet_exact)); fmpr_zero(fmprb_radref(indet_exact));
    fmpr_nan(fmprb_midref(indet_pos_rad)); fmpr_set_si(fmprb_radref(indet_pos_rad), 3);
    fmpr_nan(fmprb_midref(indet_inf_rad)); fmpr_pos_inf(fmprb_radref(indet_inf_rad));

    ASSERT(fmprb_is_zero(zero));
    ASSERT(!fmprb_is_zero(pos));
    ASSERT(!fmprb_is_zero(neg));
    ASSERT(!fmprb_is_zero(pos_inf));
    ASSERT(!fmprb_is_zero(neg_inf));
    ASSERT(!fmprb_is_zero(pos_inf_err));
    ASSERT(!fmprb_is_zero(neg_inf_err));
    ASSERT(!fmprb_is_zero(zero_pm_inf));
    ASSERT(!fmprb_is_zero(pos_pm_inf));
    ASSERT(!fmprb_is_zero(neg_pm_inf));
    ASSERT(!fmprb_is_zero(indet_exact));
    ASSERT(!fmprb_is_zero(indet_pos_rad));
    ASSERT(!fmprb_is_zero(indet_inf_rad));

    ASSERT(!fmprb_is_nonzero(zero));
    ASSERT(fmprb_is_nonzero(pos));
    ASSERT(fmprb_is_nonzero(neg));
    ASSERT(fmprb_is_nonzero(pos_inf));
    ASSERT(fmprb_is_nonzero(neg_inf));
    ASSERT(fmprb_is_nonzero(pos_inf_err));
    ASSERT(fmprb_is_nonzero(neg_inf_err));
    ASSERT(!fmprb_is_nonzero(zero_pm_inf));
    ASSERT(!fmprb_is_nonzero(pos_pm_inf));
    ASSERT(!fmprb_is_nonzero(neg_pm_inf));
    ASSERT(!fmprb_is_nonzero(indet_exact));
    ASSERT(!fmprb_is_nonzero(indet_pos_rad));
    ASSERT(!fmprb_is_nonzero(indet_inf_rad));


    ASSERT(!fmprb_is_positive(zero));
    ASSERT(fmprb_is_positive(pos));
    ASSERT(!fmprb_is_positive(neg));
    ASSERT(fmprb_is_positive(pos_inf));
    ASSERT(!fmprb_is_positive(neg_inf));
    ASSERT(fmprb_is_positive(pos_inf_err));
    ASSERT(!fmprb_is_positive(neg_inf_err));
    ASSERT(!fmprb_is_positive(zero_pm_inf));
    ASSERT(!fmprb_is_positive(pos_pm_inf));
    ASSERT(!fmprb_is_positive(neg_pm_inf));
    ASSERT(!fmprb_is_positive(indet_exact));
    ASSERT(!fmprb_is_positive(indet_pos_rad));
    ASSERT(!fmprb_is_positive(indet_inf_rad));

    ASSERT(!fmprb_is_negative(zero));
    ASSERT(!fmprb_is_negative(pos));
    ASSERT(fmprb_is_negative(neg));
    ASSERT(!fmprb_is_negative(pos_inf));
    ASSERT(fmprb_is_negative(neg_inf));
    ASSERT(!fmprb_is_negative(pos_inf_err));
    ASSERT(fmprb_is_negative(neg_inf_err));
    ASSERT(!fmprb_is_negative(zero_pm_inf));
    ASSERT(!fmprb_is_negative(pos_pm_inf));
    ASSERT(!fmprb_is_negative(neg_pm_inf));
    ASSERT(!fmprb_is_negative(indet_exact));
    ASSERT(!fmprb_is_negative(indet_pos_rad));
    ASSERT(!fmprb_is_negative(indet_inf_rad));

    ASSERT(fmprb_is_nonnegative(zero));
    ASSERT(fmprb_is_nonnegative(pos));
    ASSERT(!fmprb_is_nonnegative(neg));
    ASSERT(fmprb_is_nonnegative(pos_inf));
    ASSERT(!fmprb_is_nonnegative(neg_inf));
    ASSERT(fmprb_is_nonnegative(pos_inf_err));
    ASSERT(!fmprb_is_nonnegative(neg_inf_err));
    ASSERT(!fmprb_is_nonnegative(zero_pm_inf));
    ASSERT(!fmprb_is_nonnegative(pos_pm_inf));
    ASSERT(!fmprb_is_nonnegative(neg_pm_inf));
    ASSERT(!fmprb_is_nonnegative(indet_exact));
    ASSERT(!fmprb_is_nonnegative(indet_pos_rad));
    ASSERT(!fmprb_is_nonnegative(indet_inf_rad));

    ASSERT(fmprb_is_nonpositive(zero));
    ASSERT(!fmprb_is_nonpositive(pos));
    ASSERT(fmprb_is_nonpositive(neg));
    ASSERT(!fmprb_is_nonpositive(pos_inf));
    ASSERT(fmprb_is_nonpositive(neg_inf));
    ASSERT(!fmprb_is_nonpositive(pos_inf_err));
    ASSERT(fmprb_is_nonpositive(neg_inf_err));
    ASSERT(!fmprb_is_nonpositive(zero_pm_inf));
    ASSERT(!fmprb_is_nonpositive(pos_pm_inf));
    ASSERT(!fmprb_is_nonpositive(neg_pm_inf));
    ASSERT(!fmprb_is_nonpositive(indet_exact));
    ASSERT(!fmprb_is_nonpositive(indet_pos_rad));
    ASSERT(!fmprb_is_nonpositive(indet_inf_rad));

    ASSERT(!fmprb_contains_negative(zero));
    ASSERT(!fmprb_contains_negative(pos));
    ASSERT(fmprb_contains_negative(neg));
    ASSERT(!fmprb_contains_negative(pos_inf));
    ASSERT(fmprb_contains_negative(neg_inf));
    ASSERT(!fmprb_contains_negative(pos_inf_err));
    ASSERT(fmprb_contains_negative(neg_inf_err));
    ASSERT(fmprb_contains_negative(zero_pm_inf));
    ASSERT(fmprb_contains_negative(pos_pm_inf));
    ASSERT(fmprb_contains_negative(neg_pm_inf));
    ASSERT(fmprb_contains_negative(indet_exact));
    ASSERT(fmprb_contains_negative(indet_pos_rad));
    ASSERT(fmprb_contains_negative(indet_inf_rad));

    ASSERT(fmprb_contains_nonpositive(zero));
    ASSERT(!fmprb_contains_nonpositive(pos));
    ASSERT(fmprb_contains_nonpositive(neg));
    ASSERT(!fmprb_contains_nonpositive(pos_inf));
    ASSERT(fmprb_contains_nonpositive(neg_inf));
    ASSERT(!fmprb_contains_nonpositive(pos_inf_err));
    ASSERT(fmprb_contains_nonpositive(neg_inf_err));
    ASSERT(fmprb_contains_nonpositive(zero_pm_inf));
    ASSERT(fmprb_contains_nonpositive(pos_pm_inf));
    ASSERT(fmprb_contains_nonpositive(neg_pm_inf));
    ASSERT(fmprb_contains_nonpositive(indet_exact));
    ASSERT(fmprb_contains_nonpositive(indet_pos_rad));
    ASSERT(fmprb_contains_nonpositive(indet_inf_rad));

    ASSERT(!fmprb_contains_positive(zero));
    ASSERT(fmprb_contains_positive(pos));
    ASSERT(!fmprb_contains_positive(neg));
    ASSERT(fmprb_contains_positive(pos_inf));
    ASSERT(!fmprb_contains_positive(neg_inf));
    ASSERT(fmprb_contains_positive(pos_inf_err));
    ASSERT(!fmprb_contains_positive(neg_inf_err));
    ASSERT(fmprb_contains_positive(zero_pm_inf));
    ASSERT(fmprb_contains_positive(pos_pm_inf));
    ASSERT(fmprb_contains_positive(neg_pm_inf));
    ASSERT(fmprb_contains_positive(indet_exact));
    ASSERT(fmprb_contains_positive(indet_pos_rad));
    ASSERT(fmprb_contains_positive(indet_inf_rad));

    ASSERT(fmprb_contains_nonnegative(zero));
    ASSERT(fmprb_contains_nonnegative(pos));
    ASSERT(!fmprb_contains_nonnegative(neg));
    ASSERT(fmprb_contains_nonnegative(pos_inf));
    ASSERT(!fmprb_contains_nonnegative(neg_inf));
    ASSERT(fmprb_contains_nonnegative(pos_inf_err));
    ASSERT(!fmprb_contains_nonnegative(neg_inf_err));
    ASSERT(fmprb_contains_nonnegative(zero_pm_inf));
    ASSERT(fmprb_contains_nonnegative(pos_pm_inf));
    ASSERT(fmprb_contains_nonnegative(neg_pm_inf));
    ASSERT(fmprb_contains_nonnegative(indet_exact));
    ASSERT(fmprb_contains_nonnegative(indet_pos_rad));
    ASSERT(fmprb_contains_nonnegative(indet_inf_rad));

    ASSERT(fmprb_is_finite(zero));
    ASSERT(fmprb_is_finite(pos));
    ASSERT(fmprb_is_finite(neg));
    ASSERT(!fmprb_is_finite(pos_inf));
    ASSERT(!fmprb_is_finite(neg_inf));
    ASSERT(!fmprb_is_finite(pos_inf_err));
    ASSERT(!fmprb_is_finite(neg_inf_err));
    ASSERT(!fmprb_is_finite(zero_pm_inf));
    ASSERT(!fmprb_is_finite(pos_pm_inf));
    ASSERT(!fmprb_is_finite(neg_pm_inf));
    ASSERT(!fmprb_is_finite(indet_exact));
    ASSERT(!fmprb_is_finite(indet_pos_rad));
    ASSERT(!fmprb_is_finite(indet_inf_rad));

    ASSERT(fmprb_contains(zero, zero));
    ASSERT(!fmprb_contains(zero, pos));
    ASSERT(!fmprb_contains(zero, neg));
    ASSERT(!fmprb_contains(zero, pos_inf));
    ASSERT(!fmprb_contains(zero, neg_inf));
    ASSERT(!fmprb_contains(zero, pos_inf_err));
    ASSERT(!fmprb_contains(zero, neg_inf_err));
    ASSERT(!fmprb_contains(zero, zero_pm_inf));
    ASSERT(!fmprb_contains(zero, pos_pm_inf));
    ASSERT(!fmprb_contains(zero, neg_pm_inf));
    ASSERT(!fmprb_contains(zero, indet_exact));
    ASSERT(!fmprb_contains(zero, indet_pos_rad));
    ASSERT(!fmprb_contains(zero, indet_inf_rad));
    ASSERT(!fmprb_contains(pos, zero));
    ASSERT(fmprb_contains(pos, pos));
    ASSERT(!fmprb_contains(pos, neg));
    ASSERT(!fmprb_contains(pos, pos_inf));
    ASSERT(!fmprb_contains(pos, neg_inf));
    ASSERT(!fmprb_contains(pos, pos_inf_err));
    ASSERT(!fmprb_contains(pos, neg_inf_err));
    ASSERT(!fmprb_contains(pos, zero_pm_inf));
    ASSERT(!fmprb_contains(pos, pos_pm_inf));
    ASSERT(!fmprb_contains(pos, neg_pm_inf));
    ASSERT(!fmprb_contains(pos, indet_exact));
    ASSERT(!fmprb_contains(pos, indet_pos_rad));
    ASSERT(!fmprb_contains(pos, indet_inf_rad));
    ASSERT(!fmprb_contains(neg, zero));
    ASSERT(!fmprb_contains(neg, pos));
    ASSERT(fmprb_contains(neg, neg));
    ASSERT(!fmprb_contains(neg, pos_inf));
    ASSERT(!fmprb_contains(neg, neg_inf));
    ASSERT(!fmprb_contains(neg, pos_inf_err));
    ASSERT(!fmprb_contains(neg, neg_inf_err));
    ASSERT(!fmprb_contains(neg, zero_pm_inf));
    ASSERT(!fmprb_contains(neg, pos_pm_inf));
    ASSERT(!fmprb_contains(neg, neg_pm_inf));
    ASSERT(!fmprb_contains(neg, indet_exact));
    ASSERT(!fmprb_contains(neg, indet_pos_rad));
    ASSERT(!fmprb_contains(neg, indet_inf_rad));
    ASSERT(!fmprb_contains(pos_inf, zero));
    ASSERT(!fmprb_contains(pos_inf, pos));
    ASSERT(!fmprb_contains(pos_inf, neg));
    ASSERT(fmprb_contains(pos_inf, pos_inf));
    ASSERT(!fmprb_contains(pos_inf, neg_inf));
    ASSERT(fmprb_contains(pos_inf, pos_inf_err));
    ASSERT(!fmprb_contains(pos_inf, neg_inf_err));
    ASSERT(!fmprb_contains(pos_inf, zero_pm_inf));
    ASSERT(!fmprb_contains(pos_inf, pos_pm_inf));
    ASSERT(!fmprb_contains(pos_inf, neg_pm_inf));
    ASSERT(!fmprb_contains(pos_inf, indet_exact));
    ASSERT(!fmprb_contains(pos_inf, indet_pos_rad));
    ASSERT(!fmprb_contains(pos_inf, indet_inf_rad));
    ASSERT(!fmprb_contains(neg_inf, zero));
    ASSERT(!fmprb_contains(neg_inf, pos));
    ASSERT(!fmprb_contains(neg_inf, neg));
    ASSERT(!fmprb_contains(neg_inf, pos_inf));
    ASSERT(fmprb_contains(neg_inf, neg_inf));
    ASSERT(!fmprb_contains(neg_inf, pos_inf_err));
    ASSERT(fmprb_contains(neg_inf, neg_inf_err));
    ASSERT(!fmprb_contains(neg_inf, zero_pm_inf));
    ASSERT(!fmprb_contains(neg_inf, pos_pm_inf));
    ASSERT(!fmprb_contains(neg_inf, neg_pm_inf));
    ASSERT(!fmprb_contains(neg_inf, indet_exact));
    ASSERT(!fmprb_contains(neg_inf, indet_pos_rad));
    ASSERT(!fmprb_contains(neg_inf, indet_inf_rad));
    ASSERT(!fmprb_contains(pos_inf_err, zero));
    ASSERT(!fmprb_contains(pos_inf_err, pos));
    ASSERT(!fmprb_contains(pos_inf_err, neg));
    ASSERT(fmprb_contains(pos_inf_err, pos_inf));
    ASSERT(!fmprb_contains(pos_inf_err, neg_inf));
    ASSERT(fmprb_contains(pos_inf_err, pos_inf_err));
    ASSERT(!fmprb_contains(pos_inf_err, neg_inf_err));
    ASSERT(!fmprb_contains(pos_inf_err, zero_pm_inf));
    ASSERT(!fmprb_contains(pos_inf_err, pos_pm_inf));
    ASSERT(!fmprb_contains(pos_inf_err, neg_pm_inf));
    ASSERT(!fmprb_contains(pos_inf_err, indet_exact));
    ASSERT(!fmprb_contains(pos_inf_err, indet_pos_rad));
    ASSERT(!fmprb_contains(pos_inf_err, indet_inf_rad));
    ASSERT(!fmprb_contains(neg_inf_err, zero));
    ASSERT(!fmprb_contains(neg_inf_err, pos));
    ASSERT(!fmprb_contains(neg_inf_err, neg));
    ASSERT(!fmprb_contains(neg_inf_err, pos_inf));
    ASSERT(fmprb_contains(neg_inf_err, neg_inf));
    ASSERT(!fmprb_contains(neg_inf_err, pos_inf_err));
    ASSERT(fmprb_contains(neg_inf_err, neg_inf_err));
    ASSERT(!fmprb_contains(neg_inf_err, zero_pm_inf));
    ASSERT(!fmprb_contains(neg_inf_err, pos_pm_inf));
    ASSERT(!fmprb_contains(neg_inf_err, neg_pm_inf));
    ASSERT(!fmprb_contains(neg_inf_err, indet_exact));
    ASSERT(!fmprb_contains(neg_inf_err, indet_pos_rad));
    ASSERT(!fmprb_contains(neg_inf_err, indet_inf_rad));
    ASSERT(fmprb_contains(zero_pm_inf, zero));
    ASSERT(fmprb_contains(zero_pm_inf, pos));
    ASSERT(fmprb_contains(zero_pm_inf, neg));
    ASSERT(fmprb_contains(zero_pm_inf, pos_inf));
    ASSERT(fmprb_contains(zero_pm_inf, neg_inf));
    ASSERT(fmprb_contains(zero_pm_inf, pos_inf_err));
    ASSERT(fmprb_contains(zero_pm_inf, neg_inf_err));
    ASSERT(fmprb_contains(zero_pm_inf, zero_pm_inf));
    ASSERT(fmprb_contains(zero_pm_inf, pos_pm_inf));
    ASSERT(fmprb_contains(zero_pm_inf, neg_pm_inf));
    ASSERT(!fmprb_contains(zero_pm_inf, indet_exact));
    ASSERT(!fmprb_contains(zero_pm_inf, indet_pos_rad));
    ASSERT(!fmprb_contains(zero_pm_inf, indet_inf_rad));
    ASSERT(fmprb_contains(pos_pm_inf, zero));
    ASSERT(fmprb_contains(pos_pm_inf, pos));
    ASSERT(fmprb_contains(pos_pm_inf, neg));
    ASSERT(fmprb_contains(pos_pm_inf, pos_inf));
    ASSERT(fmprb_contains(pos_pm_inf, neg_inf));
    ASSERT(fmprb_contains(pos_pm_inf, pos_inf_err));
    ASSERT(fmprb_contains(pos_pm_inf, neg_inf_err));
    ASSERT(fmprb_contains(pos_pm_inf, zero_pm_inf));
    ASSERT(fmprb_contains(pos_pm_inf, pos_pm_inf));
    ASSERT(fmprb_contains(pos_pm_inf, neg_pm_inf));
    ASSERT(!fmprb_contains(pos_pm_inf, indet_exact));
    ASSERT(!fmprb_contains(pos_pm_inf, indet_pos_rad));
    ASSERT(!fmprb_contains(pos_pm_inf, indet_inf_rad));
    ASSERT(fmprb_contains(neg_pm_inf, zero));
    ASSERT(fmprb_contains(neg_pm_inf, pos));
    ASSERT(fmprb_contains(neg_pm_inf, neg));
    ASSERT(fmprb_contains(neg_pm_inf, pos_inf));
    ASSERT(fmprb_contains(neg_pm_inf, neg_inf));
    ASSERT(fmprb_contains(neg_pm_inf, pos_inf_err));
    ASSERT(fmprb_contains(neg_pm_inf, neg_inf_err));
    ASSERT(fmprb_contains(neg_pm_inf, zero_pm_inf));
    ASSERT(fmprb_contains(neg_pm_inf, pos_pm_inf));
    ASSERT(fmprb_contains(neg_pm_inf, neg_pm_inf));
    ASSERT(!fmprb_contains(neg_pm_inf, indet_exact));
    ASSERT(!fmprb_contains(neg_pm_inf, indet_pos_rad));
    ASSERT(!fmprb_contains(neg_pm_inf, indet_inf_rad));
    ASSERT(fmprb_contains(indet_exact, zero));
    ASSERT(fmprb_contains(indet_exact, pos));
    ASSERT(fmprb_contains(indet_exact, neg));
    ASSERT(fmprb_contains(indet_exact, pos_inf));
    ASSERT(fmprb_contains(indet_exact, neg_inf));
    ASSERT(fmprb_contains(indet_exact, pos_inf_err));
    ASSERT(fmprb_contains(indet_exact, neg_inf_err));
    ASSERT(fmprb_contains(indet_exact, zero_pm_inf));
    ASSERT(fmprb_contains(indet_exact, pos_pm_inf));
    ASSERT(fmprb_contains(indet_exact, neg_pm_inf));
    ASSERT(fmprb_contains(indet_exact, indet_exact));
    ASSERT(fmprb_contains(indet_exact, indet_pos_rad));
    ASSERT(fmprb_contains(indet_exact, indet_inf_rad));
    ASSERT(fmprb_contains(indet_pos_rad, zero));
    ASSERT(fmprb_contains(indet_pos_rad, pos));
    ASSERT(fmprb_contains(indet_pos_rad, neg));
    ASSERT(fmprb_contains(indet_pos_rad, pos_inf));
    ASSERT(fmprb_contains(indet_pos_rad, neg_inf));
    ASSERT(fmprb_contains(indet_pos_rad, pos_inf_err));
    ASSERT(fmprb_contains(indet_pos_rad, neg_inf_err));
    ASSERT(fmprb_contains(indet_pos_rad, zero_pm_inf));
    ASSERT(fmprb_contains(indet_pos_rad, pos_pm_inf));
    ASSERT(fmprb_contains(indet_pos_rad, neg_pm_inf));
    ASSERT(fmprb_contains(indet_pos_rad, indet_exact));
    ASSERT(fmprb_contains(indet_pos_rad, indet_pos_rad));
    ASSERT(fmprb_contains(indet_pos_rad, indet_inf_rad));
    ASSERT(fmprb_contains(indet_inf_rad, zero));
    ASSERT(fmprb_contains(indet_inf_rad, pos));
    ASSERT(fmprb_contains(indet_inf_rad, neg));
    ASSERT(fmprb_contains(indet_inf_rad, pos_inf));
    ASSERT(fmprb_contains(indet_inf_rad, neg_inf));
    ASSERT(fmprb_contains(indet_inf_rad, pos_inf_err));
    ASSERT(fmprb_contains(indet_inf_rad, neg_inf_err));
    ASSERT(fmprb_contains(indet_inf_rad, zero_pm_inf));
    ASSERT(fmprb_contains(indet_inf_rad, pos_pm_inf));
    ASSERT(fmprb_contains(indet_inf_rad, neg_pm_inf));
    ASSERT(fmprb_contains(indet_inf_rad, indet_exact));
    ASSERT(fmprb_contains(indet_inf_rad, indet_pos_rad));
    ASSERT(fmprb_contains(indet_inf_rad, indet_inf_rad));

    ASSERT(fmprb_overlaps(zero, zero));
    ASSERT(!fmprb_overlaps(zero, pos));
    ASSERT(!fmprb_overlaps(zero, neg));
    ASSERT(!fmprb_overlaps(zero, pos_inf));
    ASSERT(!fmprb_overlaps(zero, neg_inf));
    ASSERT(!fmprb_overlaps(zero, pos_inf_err));
    ASSERT(!fmprb_overlaps(zero, neg_inf_err));
    ASSERT(fmprb_overlaps(zero, zero_pm_inf));
    ASSERT(fmprb_overlaps(zero, pos_pm_inf));
    ASSERT(fmprb_overlaps(zero, neg_pm_inf));
    ASSERT(fmprb_overlaps(zero, indet_exact));
    ASSERT(fmprb_overlaps(zero, indet_pos_rad));
    ASSERT(fmprb_overlaps(zero, indet_inf_rad));
    ASSERT(!fmprb_overlaps(pos, zero));
    ASSERT(fmprb_overlaps(pos, pos));
    ASSERT(!fmprb_overlaps(pos, neg));
    ASSERT(!fmprb_overlaps(pos, pos_inf));
    ASSERT(!fmprb_overlaps(pos, neg_inf));
    ASSERT(!fmprb_overlaps(pos, pos_inf_err));
    ASSERT(!fmprb_overlaps(pos, neg_inf_err));
    ASSERT(fmprb_overlaps(pos, zero_pm_inf));
    ASSERT(fmprb_overlaps(pos, pos_pm_inf));
    ASSERT(fmprb_overlaps(pos, neg_pm_inf));
    ASSERT(fmprb_overlaps(pos, indet_exact));
    ASSERT(fmprb_overlaps(pos, indet_pos_rad));
    ASSERT(fmprb_overlaps(pos, indet_inf_rad));
    ASSERT(!fmprb_overlaps(neg, zero));
    ASSERT(!fmprb_overlaps(neg, pos));
    ASSERT(fmprb_overlaps(neg, neg));
    ASSERT(!fmprb_overlaps(neg, pos_inf));
    ASSERT(!fmprb_overlaps(neg, neg_inf));
    ASSERT(!fmprb_overlaps(neg, pos_inf_err));
    ASSERT(!fmprb_overlaps(neg, neg_inf_err));
    ASSERT(fmprb_overlaps(neg, zero_pm_inf));
    ASSERT(fmprb_overlaps(neg, pos_pm_inf));
    ASSERT(fmprb_overlaps(neg, neg_pm_inf));
    ASSERT(fmprb_overlaps(neg, indet_exact));
    ASSERT(fmprb_overlaps(neg, indet_pos_rad));
    ASSERT(fmprb_overlaps(neg, indet_inf_rad));
    ASSERT(!fmprb_overlaps(pos_inf, zero));
    ASSERT(!fmprb_overlaps(pos_inf, pos));
    ASSERT(!fmprb_overlaps(pos_inf, neg));
    ASSERT(fmprb_overlaps(pos_inf, pos_inf));
    ASSERT(!fmprb_overlaps(pos_inf, neg_inf));
    ASSERT(fmprb_overlaps(pos_inf, pos_inf_err));
    ASSERT(!fmprb_overlaps(pos_inf, neg_inf_err));
    ASSERT(fmprb_overlaps(pos_inf, zero_pm_inf));
    ASSERT(fmprb_overlaps(pos_inf, pos_pm_inf));
    ASSERT(fmprb_overlaps(pos_inf, neg_pm_inf));
    ASSERT(fmprb_overlaps(pos_inf, indet_exact));
    ASSERT(fmprb_overlaps(pos_inf, indet_pos_rad));
    ASSERT(fmprb_overlaps(pos_inf, indet_inf_rad));
    ASSERT(!fmprb_overlaps(neg_inf, zero));
    ASSERT(!fmprb_overlaps(neg_inf, pos));
    ASSERT(!fmprb_overlaps(neg_inf, neg));
    ASSERT(!fmprb_overlaps(neg_inf, pos_inf));
    ASSERT(fmprb_overlaps(neg_inf, neg_inf));
    ASSERT(!fmprb_overlaps(neg_inf, pos_inf_err));
    ASSERT(fmprb_overlaps(neg_inf, neg_inf_err));
    ASSERT(fmprb_overlaps(neg_inf, zero_pm_inf));
    ASSERT(fmprb_overlaps(neg_inf, pos_pm_inf));
    ASSERT(fmprb_overlaps(neg_inf, neg_pm_inf));
    ASSERT(fmprb_overlaps(neg_inf, indet_exact));
    ASSERT(fmprb_overlaps(neg_inf, indet_pos_rad));
    ASSERT(fmprb_overlaps(neg_inf, indet_inf_rad));
    ASSERT(!fmprb_overlaps(pos_inf_err, zero));
    ASSERT(!fmprb_overlaps(pos_inf_err, pos));
    ASSERT(!fmprb_overlaps(pos_inf_err, neg));
    ASSERT(fmprb_overlaps(pos_inf_err, pos_inf));
    ASSERT(!fmprb_overlaps(pos_inf_err, neg_inf));
    ASSERT(fmprb_overlaps(pos_inf_err, pos_inf_err));
    ASSERT(!fmprb_overlaps(pos_inf_err, neg_inf_err));
    ASSERT(fmprb_overlaps(pos_inf_err, zero_pm_inf));
    ASSERT(fmprb_overlaps(pos_inf_err, pos_pm_inf));
    ASSERT(fmprb_overlaps(pos_inf_err, neg_pm_inf));
    ASSERT(fmprb_overlaps(pos_inf_err, indet_exact));
    ASSERT(fmprb_overlaps(pos_inf_err, indet_pos_rad));
    ASSERT(fmprb_overlaps(pos_inf_err, indet_inf_rad));
    ASSERT(!fmprb_overlaps(neg_inf_err, zero));
    ASSERT(!fmprb_overlaps(neg_inf_err, pos));
    ASSERT(!fmprb_overlaps(neg_inf_err, neg));
    ASSERT(!fmprb_overlaps(neg_inf_err, pos_inf));
    ASSERT(fmprb_overlaps(neg_inf_err, neg_inf));
    ASSERT(!fmprb_overlaps(neg_inf_err, pos_inf_err));
    ASSERT(fmprb_overlaps(neg_inf_err, neg_inf_err));
    ASSERT(fmprb_overlaps(neg_inf_err, zero_pm_inf));
    ASSERT(fmprb_overlaps(neg_inf_err, pos_pm_inf));
    ASSERT(fmprb_overlaps(neg_inf_err, neg_pm_inf));
    ASSERT(fmprb_overlaps(neg_inf_err, indet_exact));
    ASSERT(fmprb_overlaps(neg_inf_err, indet_pos_rad));
    ASSERT(fmprb_overlaps(neg_inf_err, indet_inf_rad));
    ASSERT(fmprb_overlaps(zero_pm_inf, zero));
    ASSERT(fmprb_overlaps(zero_pm_inf, pos));
    ASSERT(fmprb_overlaps(zero_pm_inf, neg));
    ASSERT(fmprb_overlaps(zero_pm_inf, pos_inf));
    ASSERT(fmprb_overlaps(zero_pm_inf, neg_inf));
    ASSERT(fmprb_overlaps(zero_pm_inf, pos_inf_err));
    ASSERT(fmprb_overlaps(zero_pm_inf, neg_inf_err));
    ASSERT(fmprb_overlaps(zero_pm_inf, zero_pm_inf));
    ASSERT(fmprb_overlaps(zero_pm_inf, pos_pm_inf));
    ASSERT(fmprb_overlaps(zero_pm_inf, neg_pm_inf));
    ASSERT(fmprb_overlaps(zero_pm_inf, indet_exact));
    ASSERT(fmprb_overlaps(zero_pm_inf, indet_pos_rad));
    ASSERT(fmprb_overlaps(zero_pm_inf, indet_inf_rad));
    ASSERT(fmprb_overlaps(pos_pm_inf, zero));
    ASSERT(fmprb_overlaps(pos_pm_inf, pos));
    ASSERT(fmprb_overlaps(pos_pm_inf, neg));
    ASSERT(fmprb_overlaps(pos_pm_inf, pos_inf));
    ASSERT(fmprb_overlaps(pos_pm_inf, neg_inf));
    ASSERT(fmprb_overlaps(pos_pm_inf, pos_inf_err));
    ASSERT(fmprb_overlaps(pos_pm_inf, neg_inf_err));
    ASSERT(fmprb_overlaps(pos_pm_inf, zero_pm_inf));
    ASSERT(fmprb_overlaps(pos_pm_inf, pos_pm_inf));
    ASSERT(fmprb_overlaps(pos_pm_inf, neg_pm_inf));
    ASSERT(fmprb_overlaps(pos_pm_inf, indet_exact));
    ASSERT(fmprb_overlaps(pos_pm_inf, indet_pos_rad));
    ASSERT(fmprb_overlaps(pos_pm_inf, indet_inf_rad));
    ASSERT(fmprb_overlaps(neg_pm_inf, zero));
    ASSERT(fmprb_overlaps(neg_pm_inf, pos));
    ASSERT(fmprb_overlaps(neg_pm_inf, neg));
    ASSERT(fmprb_overlaps(neg_pm_inf, pos_inf));
    ASSERT(fmprb_overlaps(neg_pm_inf, neg_inf));
    ASSERT(fmprb_overlaps(neg_pm_inf, pos_inf_err));
    ASSERT(fmprb_overlaps(neg_pm_inf, neg_inf_err));
    ASSERT(fmprb_overlaps(neg_pm_inf, zero_pm_inf));
    ASSERT(fmprb_overlaps(neg_pm_inf, pos_pm_inf));
    ASSERT(fmprb_overlaps(neg_pm_inf, neg_pm_inf));
    ASSERT(fmprb_overlaps(neg_pm_inf, indet_exact));
    ASSERT(fmprb_overlaps(neg_pm_inf, indet_pos_rad));
    ASSERT(fmprb_overlaps(neg_pm_inf, indet_inf_rad));
    ASSERT(fmprb_overlaps(indet_exact, zero));
    ASSERT(fmprb_overlaps(indet_exact, pos));
    ASSERT(fmprb_overlaps(indet_exact, neg));
    ASSERT(fmprb_overlaps(indet_exact, pos_inf));
    ASSERT(fmprb_overlaps(indet_exact, neg_inf));
    ASSERT(fmprb_overlaps(indet_exact, pos_inf_err));
    ASSERT(fmprb_overlaps(indet_exact, neg_inf_err));
    ASSERT(fmprb_overlaps(indet_exact, zero_pm_inf));
    ASSERT(fmprb_overlaps(indet_exact, pos_pm_inf));
    ASSERT(fmprb_overlaps(indet_exact, neg_pm_inf));
    ASSERT(fmprb_overlaps(indet_exact, indet_exact));
    ASSERT(fmprb_overlaps(indet_exact, indet_pos_rad));
    ASSERT(fmprb_overlaps(indet_exact, indet_inf_rad));
    ASSERT(fmprb_overlaps(indet_pos_rad, zero));
    ASSERT(fmprb_overlaps(indet_pos_rad, pos));
    ASSERT(fmprb_overlaps(indet_pos_rad, neg));
    ASSERT(fmprb_overlaps(indet_pos_rad, pos_inf));
    ASSERT(fmprb_overlaps(indet_pos_rad, neg_inf));
    ASSERT(fmprb_overlaps(indet_pos_rad, pos_inf_err));
    ASSERT(fmprb_overlaps(indet_pos_rad, neg_inf_err));
    ASSERT(fmprb_overlaps(indet_pos_rad, zero_pm_inf));
    ASSERT(fmprb_overlaps(indet_pos_rad, pos_pm_inf));
    ASSERT(fmprb_overlaps(indet_pos_rad, neg_pm_inf));
    ASSERT(fmprb_overlaps(indet_pos_rad, indet_exact));
    ASSERT(fmprb_overlaps(indet_pos_rad, indet_pos_rad));
    ASSERT(fmprb_overlaps(indet_pos_rad, indet_inf_rad));
    ASSERT(fmprb_overlaps(indet_inf_rad, zero));
    ASSERT(fmprb_overlaps(indet_inf_rad, pos));
    ASSERT(fmprb_overlaps(indet_inf_rad, neg));
    ASSERT(fmprb_overlaps(indet_inf_rad, pos_inf));
    ASSERT(fmprb_overlaps(indet_inf_rad, neg_inf));
    ASSERT(fmprb_overlaps(indet_inf_rad, pos_inf_err));
    ASSERT(fmprb_overlaps(indet_inf_rad, neg_inf_err));
    ASSERT(fmprb_overlaps(indet_inf_rad, zero_pm_inf));
    ASSERT(fmprb_overlaps(indet_inf_rad, pos_pm_inf));
    ASSERT(fmprb_overlaps(indet_inf_rad, neg_pm_inf));
    ASSERT(fmprb_overlaps(indet_inf_rad, indet_exact));
    ASSERT(fmprb_overlaps(indet_inf_rad, indet_pos_rad));
    ASSERT(fmprb_overlaps(indet_inf_rad, indet_inf_rad));

    {
        fmpz_t t;
        fmpz_init(t);

        ASSERT(fmprb_get_unique_fmpz(t, zero));
        ASSERT(!fmprb_get_unique_fmpz(t, pos_inf));
        ASSERT(!fmprb_get_unique_fmpz(t, neg_inf));
        ASSERT(!fmprb_get_unique_fmpz(t, pos_inf_err));
        ASSERT(!fmprb_get_unique_fmpz(t, neg_inf_err));
        ASSERT(!fmprb_get_unique_fmpz(t, zero_pm_inf));
        ASSERT(!fmprb_get_unique_fmpz(t, pos_pm_inf));
        ASSERT(!fmprb_get_unique_fmpz(t, neg_pm_inf));
        ASSERT(!fmprb_get_unique_fmpz(t, indet_exact));
        ASSERT(!fmprb_get_unique_fmpz(t, indet_pos_rad));
        ASSERT(!fmprb_get_unique_fmpz(t, indet_inf_rad));

        fmpz_clear(t);
    }

    {
        fmpr_t b;
        long wp = 30;

        fmpr_init(b);

        fmprb_get_abs_ubound_fmpr(b, zero, wp); ASSERT(fmpr_is_zero(b));
        fmprb_get_abs_ubound_fmpr(b, pos, wp); ASSERT(fmpr_sgn(b) > 0);
        fmprb_get_abs_ubound_fmpr(b, neg, wp); ASSERT(fmpr_sgn(b) > 0);
        fmprb_get_abs_ubound_fmpr(b, pos_inf, wp); ASSERT(fmpr_is_pos_inf(b));
        fmprb_get_abs_ubound_fmpr(b, neg_inf, wp); ASSERT(fmpr_is_pos_inf(b));
        fmprb_get_abs_ubound_fmpr(b, pos_inf_err, wp); ASSERT(fmpr_is_pos_inf(b));
        fmprb_get_abs_ubound_fmpr(b, neg_inf_err, wp); ASSERT(fmpr_is_pos_inf(b));
        fmprb_get_abs_ubound_fmpr(b, zero_pm_inf, wp); ASSERT(fmpr_is_pos_inf(b));
        fmprb_get_abs_ubound_fmpr(b, indet_exact, wp); ASSERT(fmpr_is_nan(b));
        fmprb_get_abs_ubound_fmpr(b, indet_pos_rad, wp); ASSERT(fmpr_is_nan(b));
        fmprb_get_abs_ubound_fmpr(b, indet_inf_rad, wp); ASSERT(fmpr_is_nan(b));

        fmprb_get_abs_lbound_fmpr(b, zero, wp); ASSERT(fmpr_is_zero(b));
        fmprb_get_abs_lbound_fmpr(b, pos, wp); ASSERT(fmpr_sgn(b) > 0);
        fmprb_get_abs_lbound_fmpr(b, neg, wp); ASSERT(fmpr_sgn(b) > 0);
        fmprb_get_abs_lbound_fmpr(b, pos_inf, wp); ASSERT(fmpr_is_pos_inf(b));
        fmprb_get_abs_lbound_fmpr(b, neg_inf, wp); ASSERT(fmpr_is_pos_inf(b));
        fmprb_get_abs_lbound_fmpr(b, pos_inf_err, wp); ASSERT(fmpr_is_pos_inf(b));
        fmprb_get_abs_lbound_fmpr(b, neg_inf_err, wp); ASSERT(fmpr_is_pos_inf(b));
        fmprb_get_abs_lbound_fmpr(b, zero_pm_inf, wp); ASSERT(fmpr_is_zero(b));
        fmprb_get_abs_lbound_fmpr(b, indet_exact, wp); ASSERT(fmpr_is_nan(b));
        fmprb_get_abs_lbound_fmpr(b, indet_pos_rad, wp); ASSERT(fmpr_is_nan(b));
        fmprb_get_abs_lbound_fmpr(b, indet_inf_rad, wp); ASSERT(fmpr_is_nan(b));

        fmpr_clear(b);
    }

    fmprb_clear(zero);
    fmprb_clear(pos);
    fmprb_clear(neg);
    fmprb_clear(pos_inf);
    fmprb_clear(neg_inf);
    fmprb_clear(pos_inf_err);
    fmprb_clear(neg_inf_err);
    fmprb_clear(zero_pm_inf);
    fmprb_clear(pos_pm_inf);
    fmprb_clear(neg_pm_inf);
    fmprb_clear(indet_exact);
    fmprb_clear(indet_pos_rad);
    fmprb_clear(indet_inf_rad);

    flint_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}

