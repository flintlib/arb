/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_calc.h"

/* Absolute value function on R extended to a holomorphic function in the left
   and right half planes. */
void
acb_holomorphic_abs(acb_ptr res, const acb_t z, int holomorphic, slong prec)
{
    if (!acb_is_finite(z) || (holomorphic && arb_contains_zero(acb_realref(z))))
    {
        acb_indeterminate(res);
    }
    else
    {
        if (arb_is_nonnegative(acb_realref(z)))
        {
            acb_set_round(res, z, prec);
        }
        else if (arb_is_negative(acb_realref(z)))
        {
            acb_neg_round(res, z, prec);
        }
        else
        {
            acb_t t;
            acb_init(t);
            acb_neg(t, res);
            acb_union(res, z, t, prec);
            acb_clear(t);
        }
    }
}

/* Square root function on C with detection of the branch cut. */
void
acb_holomorphic_sqrt(acb_ptr res, const acb_t z, int holomorphic, slong prec)
{
    if (!acb_is_finite(z) || (holomorphic &&
        arb_contains_zero(acb_imagref(z)) &&
        arb_contains_nonpositive(acb_realref(z))))
    {
        acb_indeterminate(res);
    }
    else
    {
        acb_sqrt(res, z, prec);
    }
}

int
f_zero(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_zero(res);

    return 0;
}

int
f_indet(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_indeterminate(res);

    return 0;
}


int
f_cube(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_cube(res, z, prec);

    return 0;
}

int
f_sin(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_sin(res, z, prec);

    return 0;
}

int
f_circle(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_one(res);
    acb_submul(res, z, z, prec);
    acb_holomorphic_sqrt(res, res, order != 0, prec);

    /* Rounding could give |z| = 1 + eps near the endpoints, but we assume
       that the interval is [-1,1] which really makes f real.  */
    if (order == 0)
        arb_zero(acb_imagref(res));

    return 0;
}

/* f(z) = sin(z + exp(z)) -- Rump's oscillatory example */
int
f_rump(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_exp(res, z, prec);
    acb_add(res, res, z, prec);
    acb_sin(res, res, prec);

    return 0;
}

/* f(z) = |z^4 + 10z^3 + 19z^2 - 6z - 6| exp(z)   (for real z) --
   Helfgott's integral on MathOverflow */
int
f_helfgott(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_add_si(res, z, 10, prec);
    acb_mul(res, res, z, prec);
    acb_add_si(res, res, 19, prec);
    acb_mul(res, res, z, prec);
    acb_add_si(res, res, -6, prec);
    acb_mul(res, res, z, prec);
    acb_add_si(res, res, -6, prec);

    acb_holomorphic_abs(res, res, order != 0, prec);

    if (acb_is_finite(res))
    {
        acb_t t;
        acb_init(t);
        acb_exp(t, z, prec);
        acb_mul(res, res, t, prec);
        acb_clear(t);
    }

    return 0;
}

/* f(z) = z sin(1/z), assume on real interval */
int
f_essing2(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    if ((order == 0) && acb_is_real(z) && arb_contains_zero(acb_realref(z)))
    {
        /* todo: arb_zero_pm_one, arb_unit_interval? */
        acb_zero(res);
        mag_one(arb_radref(acb_realref(res)));
    }
    else
    {
        acb_inv(res, z, prec);
        acb_sin(res, res, prec);
    }

    acb_mul(res, res, z, prec);

    return 0;
}

/* f(z) = sech(10(x-0.2))^2 + sech(100(x-0.4))^4 + sech(1000(x-0.6))^6 */
int
f_wolfram(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t a, b, c;

    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_init(a);
    acb_init(b);
    acb_init(c);

    acb_mul_ui(a, z, 10, prec);
    acb_sub_ui(a, a, 2, prec);
    acb_sech(a, a, prec);
    acb_pow_ui(a, a, 2, prec);

    acb_mul_ui(b, z, 100, prec);
    acb_sub_ui(b, b, 40, prec);
    acb_sech(b, b, prec);
    acb_pow_ui(b, b, 4, prec);

    acb_mul_ui(c, z, 1000, prec);
    acb_sub_ui(c, c, 600, prec);
    acb_sech(c, c, prec);
    acb_pow_ui(c, c, 6, prec);

    acb_add(res, a, b, prec);
    acb_add(res, res, c, prec);

    acb_clear(a);
    acb_clear(b);
    acb_clear(c);

    return 0;
}

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("integrate....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000 * arb_test_multiplier(); iter++)
    {
        acb_t ans, res, a, b, t;
        slong goal, prec;
        mag_t tol;
        acb_calc_integrate_opt_t opt;
        int integral;

        acb_init(ans);
        acb_init(res);
        acb_init(a);
        acb_init(b);
        acb_init(t);
        mag_init(tol);
        acb_calc_integrate_opt_init(opt);

        goal = 2 + n_randint(state, 300);
        prec = 2 + n_randint(state, 300);
        mag_set_ui_2exp_si(tol, n_randint(state, 2), -n_randint(state, 300));

        opt->eval_limit = n_randint(state, 10000);

        if (n_randint(state, 2))
            opt->depth_limit = n_randint(state, 100);

        if (n_randint(state, 2))
            opt->deg_limit = n_randint(state, 100);

        opt->use_heap = n_randint(state, 2);

        integral = n_randint(state, 9);

        if (integral == 0)
        {
            acb_randtest(a, state, 1 + n_randint(state, 200), 2);
            acb_randtest(b, state, 1 + n_randint(state, 200), 2);
            acb_cos(ans, a, prec);
            acb_cos(res, b, prec);
            acb_sub(ans, ans, res, prec);
            acb_calc_integrate(res, f_sin, NULL, a, b, goal, tol, opt, prec);
        }
        else if (integral == 1)
        {
            acb_randtest(a, state, 1 + n_randint(state, 200), 2);
            acb_randtest(b, state, 1 + n_randint(state, 200), 2);
            acb_pow_ui(ans, a, 4, prec);
            acb_pow_ui(res, b, 4, prec);
            acb_sub(ans, res, ans, prec);
            acb_mul_2exp_si(ans, ans, -2);
            acb_calc_integrate(res, f_cube, NULL, a, b, goal, tol, opt, prec);
        }
        else if (integral == 2)
        {
            acb_zero(a);
            acb_one(b);
            acb_const_pi(ans, prec);
            acb_mul_2exp_si(ans, ans, -2);
            if (n_randint(state, 2))
            {
                acb_swap(a, b);
                acb_neg(ans, ans);
            }
            acb_calc_integrate(res, f_circle, NULL, a, b, goal, tol, opt, prec);
        }
        else if (integral == 3)
        {
            acb_zero(a);
            acb_one(b);
            acb_zero(ans);
            arb_set_str(acb_realref(ans), "0.3785300171241613098817352756283519095343133 +/- 1e-40", prec);
            if (n_randint(state, 2))
            {
                acb_swap(a, b);
                acb_neg(ans, ans);
            }
            acb_calc_integrate(res, f_essing2, NULL, a, b, goal, tol, opt, prec);
        }
        else if (integral == 4)
        {
            acb_zero(a);
            acb_one(b);
            acb_zero(ans);
            arb_set_str(acb_realref(ans), "11.147310550057139733915902084255301415775813549800589418261584268232061665808482234384871404010464 +/- 2.28e-97", prec);
            if (n_randint(state, 2))
            {
                acb_swap(a, b);
                acb_neg(ans, ans);
            }
            acb_calc_integrate(res, f_helfgott, NULL, a, b, goal, tol, opt, prec);
        }
        else if (integral == 5)
        {
            acb_zero(a);
            acb_set_ui(b, 8);
            acb_zero(ans);
            arb_set_str(acb_realref(ans), "0.34740017265724780787951215911989312465745625486618018388549271361674821398878532052968510434660 +/- 5.97e-96", prec);
            if (n_randint(state, 2))
            {
                acb_swap(a, b);
                acb_neg(ans, ans);
            }
            acb_calc_integrate(res, f_rump, NULL, a, b, goal, tol, opt, prec);
        }
        else if (integral == 6)
        {
            acb_zero(a);
            acb_one(b);
            acb_zero(ans);
            arb_set_str(acb_realref(ans), "0.21080273550054927737564325570572915436090918643678119034785050587872061312814550020505868926155764 +/- 3.72e-99", prec);
            if (n_randint(state, 2))
            {
                acb_swap(a, b);
                acb_neg(ans, ans);
            }
            acb_calc_integrate(res, f_wolfram, NULL, a, b, goal, tol, opt, prec);
        }
        else if (integral == 7)
        {
            acb_randtest_special(a, state, 1 + n_randint(state, 200), 2);
            acb_randtest_special(b, state, 1 + n_randint(state, 200), 2);
            acb_zero(ans);
            acb_calc_integrate(res, f_zero, NULL, a, b, goal, tol, opt, prec);
        }
        else
        {
            acb_randtest_special(a, state, 1 + n_randint(state, 200), 2);
            acb_randtest_special(b, state, 1 + n_randint(state, 200), 2);
            acb_indeterminate(ans);
            acb_calc_integrate(res, f_indet, NULL, a, b, goal, tol, opt, prec);
        }

        if (!acb_overlaps(res, ans))
        {
            flint_printf("FAIL! (iter = %wd)\n", iter);
            flint_printf("integral = %d, prec = %wd, goal = %wd\n", integral, prec, goal);
            flint_printf("a = "); acb_printn(a, 150, 0); flint_printf("\n");
            flint_printf("b = "); acb_printn(b, 150, 0); flint_printf("\n");
            flint_printf("res = "); acb_printn(res, 150, 0); flint_printf("\n\n");
            flint_printf("ans = "); acb_printn(ans, 150, 0); flint_printf("\n\n");
            flint_abort();
        }

        acb_clear(ans);
        acb_clear(res);
        acb_clear(a);
        acb_clear(b);
        acb_clear(t);
        mag_clear(tol);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

