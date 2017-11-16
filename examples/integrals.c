/*
  Rigorous numerical integration (with fast convergence for piecewise
  holomorphic functions) using Gauss-Legendre quadrature and adaptive
  subdivision.

  Author: Fredrik Johansson.
  This file is in the public domain.
*/

#include <string.h>
#include "flint/profiler.h"
#include "arb_hypgeom.h"
#include "acb_hypgeom.h"
#include "acb_dirichlet.h"
#include "acb_calc.h"

/* ------------------------------------------------------------------------- */
/*  Useful helper functions                                                  */
/* ------------------------------------------------------------------------- */

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

/* Floor function on R extended to a piecewise holomorphic function in
   vertical strips. */
void
acb_holomorphic_floor(acb_ptr res, const acb_t z, int holomorphic, slong prec)
{
    if (!acb_is_finite(z) || (holomorphic && arb_contains_int(acb_realref(z))))
    {
        acb_indeterminate(res);
    }
    else
    {
        arb_floor(acb_realref(res), acb_realref(z), prec);
        arb_set_round(acb_imagref(res), acb_imagref(z), prec);
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

/* ------------------------------------------------------------------------- */
/*  Example integrands                                                       */
/* ------------------------------------------------------------------------- */

/* f(z) = sin(z) */
int
f_sin(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_sin(res, z, prec);

    return 0;
}

/* f(z) = floor(z) */
int
f_floor(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_holomorphic_floor(res, z, order != 0, prec);

    return 0;
}

/* f(z) = sqrt(1-z^2) */
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

/* f(z) = 1/(1+z^2) */
int
f_atanderiv(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_mul(res, z, z, prec);
    acb_add_ui(res, res, 1, prec);
    acb_inv(res, res, prec);

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

/* f(z) = |z^4 + 10z^3 + 19z^2 + 6z - 6| exp(z)   (for real z) --
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

/* f(z) = zeta(z) */
int
f_zeta(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_zeta(res, z, prec);

    return 0;
}

/* f(z) = z sin(1/z), assume on real interval */
int
f_essing(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
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

/* f(z) = exp(-z) z^1000 */
int
f_factorial1000(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t t;

    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_init(t);
    acb_pow_ui(t, z, 1000, prec);
    acb_neg(res, z);
    acb_exp(res, res, prec);
    acb_mul(res, res, t, prec);
    acb_clear(t);

    return 0;
}

/* f(z) = gamma(z) */
int
f_gamma(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_gamma(res, z, prec);

    return 0;
}

/* f(z) = sin(z) + exp(-200-z^2) */
int
f_sin_plus_small(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t t;

    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_init(t);
    acb_mul(t, z, z, prec);
    acb_add_ui(t, t, 200, prec);
    acb_neg(t, t);
    acb_exp(t, t, prec);
    acb_sin(res, z, prec);
    acb_add(res, res, t, prec);
    acb_clear(t);

    return 0;
}

/* f(z) = exp(z) */
int
f_exp(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_exp(res, z, prec);

    return 0;
}

/* f(z) = exp(-z^2) */
int
f_gaussian(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_mul(res, z, z, prec);
    acb_neg(res, res);
    acb_exp(res, res, prec);

    return 0;
}

int
f_monster(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    acb_t t;

    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_init(t);

    acb_exp(t, z, prec);
    acb_holomorphic_floor(res, t, order != 0, prec);

    if (acb_is_finite(res))
    {
        acb_sub(res, t, res, prec);
        acb_add(t, t, z, prec);
        acb_sin(t, t, prec);
        acb_mul(res, res, t, prec);
    }

    acb_clear(t);

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

/* f(z) = sech(z) */
int
f_sech(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_sech(res, z, prec);

    return 0;
}

/* f(z) = sech^3(z) */
int
f_sech3(acb_ptr res, const acb_t z, void * param, slong order, slong prec)
{
    if (order > 1)
        flint_abort();  /* Would be needed for Taylor method. */

    acb_sech(res, z, prec);
    acb_cube(res, res, prec);

    return 0;
}

/* ------------------------------------------------------------------------- */
/*  Main test program                                                        */
/* ------------------------------------------------------------------------- */

#define NUM_INTEGRALS 18

const char * descr[NUM_INTEGRALS] =
{
    "int_0^100 sin(x) dx",
    "4 int_0^1 1/(1+x^2) dx",
    "2 int_0^{inf} 1/(1+x^2) dx   (using domain truncation)",
    "4 int_0^1 sqrt(1-x^2) dx",
    "int_0^8 sin(x+exp(x)) dx",
    "int_0^100 floor(x) dx",
    "int_0^1 |x^4+10x^3+19x^2-6x-6| exp(x) dx",
    "1/(2 pi i) int zeta(s) ds  (closed path around s = 1)",
    "int_0^1 x sin(1/x) dx  (slow convergence, use -heap and/or -tol)",
    "int_0^10000 x^1000 exp(-x) dx",
    "int_1^{1+1000i} gamma(x) dx",
    "int_{-10}^{10} sin(x) + exp(-200-x^2) dx",
    "int_{-1020}^{-1000} exp(x) dx  (use -tol 0 for relative error)",
    "int_0^{inf} exp(-x^2) dx   (using domain truncation)",
    "int_0^1 sech(10(x-0.2))^2 + sech(100(x-0.4))^4 + sech(1000(x-0.6))^6 dx",
    "int_0^8 (exp(x)-floor(exp(x))) sin(x+exp(x)) dx  (use higher -eval)",
    "int_0^{inf} sech(x) dx   (using domain truncation)",
    "int_0^{inf} sech^3(x) dx   (using domain truncation)",
};

int main(int argc, char *argv[])
{
    acb_t s, t, a, b;
    mag_t tol;
    slong prec, goal, deg_limit, eval_limit, depth_limit;
    int integral, ifrom, ito;
    int i, twice, havetol, flags;

    flint_printf("Compute integrals using subdivision and Gauss-Legendre quadrature.\n");
    flint_printf("Usage: quadrature [-i n] [-prec p] [-tol eps] [-twice]\n\n");
    flint_printf("-i n       - compute integral In (0 <= n <= %d)\n", NUM_INTEGRALS - 1);
    flint_printf("-prec p    - precision in bits (default p = 333)\n");
    flint_printf("-tol eps   - approximate absolute error goal (default 2^-p)\n");
    flint_printf("-twice     - run twice (to see overhead of computing nodes)\n");
    flint_printf("\n\n");

    prec = 333;
    twice = 0;
    ifrom = 0;
    ito = NUM_INTEGRALS - 1;
    havetol = 0;
    deg_limit = -1;
    eval_limit = -1;
    depth_limit = -1;
    flags = 0;

    acb_init(a);
    acb_init(b);
    acb_init(s);
    acb_init(t);
    mag_init(tol);

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-prec"))
        {
            prec = atol(argv[i+1]);
        }
        else if (!strcmp(argv[i], "-twice"))
        {
            twice = 1;
        }
        else if (!strcmp(argv[i], "-tol"))
        {
            arb_t x;
            arb_init(x);
            arb_set_str(x, argv[i+1], 10);
            arb_get_mag(tol, x);
            arb_clear(x);
            havetol = 1;
        }
        else if (!strcmp(argv[i], "-i"))
        {
            ifrom = ito = atol(argv[i+1]);
            if (ito < 0 || ito >= NUM_INTEGRALS)
                flint_abort();
        }
        else if (!strcmp(argv[i], "-deg"))
        {
            deg_limit = atol(argv[i+1]);
        }
        else if (!strcmp(argv[i], "-eval"))
        {
            eval_limit = atol(argv[i+1]);
        }
        else if (!strcmp(argv[i], "-depth"))
        {
            depth_limit = atol(argv[i+1]);
        }
        else if (!strcmp(argv[i], "-verbose"))
        {
            flags |= ACB_CALC_VERBOSE;
        }
        else if (!strcmp(argv[i], "-verbose2"))
        {
            flags |= ACB_CALC_VERY_VERBOSE;
        }
        else if (!strcmp(argv[i], "-heap"))
        {
            flags |= ACB_CALC_INTEGRATE_HEAP;
        }
    }

    goal = prec;

    if (!havetol)
        mag_set_ui_2exp_si(tol, 1, -prec);

    for (integral = ifrom; integral <= ito; integral++)
    {
        flint_printf("I%d = %s ...\n", integral, descr[integral]);

        for (i = 0; i < 1 + twice; i++)
        {
            TIMEIT_ONCE_START
            switch (integral)
            {
            case 0:
                acb_set_d(a, 0);
                acb_set_d(b, 100);
                acb_calc_integrate(s, f_sin, NULL, a, b, goal, tol, deg_limit, eval_limit, depth_limit, flags, prec);
                break;

            case 1:
                acb_set_d(a, 0);
                acb_set_d(b, 1);
                acb_calc_integrate(s, f_atanderiv, NULL, a, b, goal, tol, deg_limit, eval_limit, depth_limit, flags, prec);
                acb_mul_2exp_si(s, s, 2);
                break;

            case 2:
                acb_set_d(a, 0);
                acb_one(b);
                acb_mul_2exp_si(b, b, goal);
                acb_calc_integrate(s, f_atanderiv, NULL, a, b, goal, tol, deg_limit, eval_limit, depth_limit, flags, prec);
                arb_add_error_2exp_si(acb_realref(s), -goal);
                acb_mul_2exp_si(s, s, 1);
                break;

            case 3:
                acb_set_d(a, 0);
                acb_set_d(b, 1);
                acb_calc_integrate(s, f_circle, NULL, a, b, goal, tol, deg_limit, eval_limit, depth_limit, flags, prec);
                acb_mul_2exp_si(s, s, 2);
                break;

            case 4:
                acb_set_d(a, 0);
                acb_set_d(b, 8);
                acb_calc_integrate(s, f_rump, NULL, a, b, goal, tol, deg_limit, eval_limit, depth_limit, flags, prec);
                break;

            case 5:
                acb_set_d(a, 0);
                acb_set_d(b, 100);
                acb_calc_integrate(s, f_floor, NULL, a, b, goal, tol, deg_limit, eval_limit, depth_limit, flags, prec);
                break;

            case 6:
                acb_set_d(a, 0);
                acb_set_d(b, 1);
                acb_calc_integrate(s, f_helfgott, NULL, a, b, goal, tol, deg_limit, eval_limit, depth_limit, flags, prec);
                break;

            case 7:
                acb_zero(s);

                acb_set_d_d(a, -1.0, -1.0);
                acb_set_d_d(b, 2.0, -1.0);
                acb_calc_integrate(t, f_zeta, NULL, a, b, goal, tol, deg_limit, eval_limit, depth_limit, flags, prec);
                acb_add(s, s, t, prec);

                acb_set_d_d(a, 2.0, -1.0);
                acb_set_d_d(b, 2.0, 1.0);
                acb_calc_integrate(t, f_zeta, NULL, a, b, goal, tol, deg_limit, eval_limit, depth_limit, flags, prec);
                acb_add(s, s, t, prec);

                acb_set_d_d(a, 2.0, 1.0);
                acb_set_d_d(b, -1.0, 1.0);
                acb_calc_integrate(t, f_zeta, NULL, a, b, goal, tol, deg_limit, eval_limit, depth_limit, flags, prec);
                acb_add(s, s, t, prec);

                acb_set_d_d(a, -1.0, 1.0);
                acb_set_d_d(b, -1.0, -1.0);
                acb_calc_integrate(t, f_zeta, NULL, a, b, goal, tol, deg_limit, eval_limit, depth_limit, flags, prec);
                acb_add(s, s, t, prec);

                acb_const_pi(t, prec);
                acb_div(s, s, t, prec);
                acb_mul_2exp_si(s, s, -1);
                acb_div_onei(s, s);
                break;

            case 8:
                acb_set_d(a, 0);
                acb_set_d(b, 1);
                acb_calc_integrate(s, f_essing, NULL, a, b, goal, tol, deg_limit, eval_limit, depth_limit, flags, prec);
                break;

            case 9:
                acb_set_d(a, 0);
                acb_set_d(b, 10000);
                acb_calc_integrate(s, f_factorial1000, NULL, a, b, goal, tol, deg_limit, eval_limit, depth_limit, flags, prec);
                break;

            case 10:
                acb_set_d_d(a, 1.0, 0.0);
                acb_set_d_d(b, 1.0, 1000.0);
                acb_calc_integrate(s, f_gamma, NULL, a, b, goal, tol, deg_limit, eval_limit, depth_limit, flags, prec);
                break;

            case 11:
                acb_set_d(a, -10.0);
                acb_set_d(b, 10.0);
                acb_calc_integrate(s, f_sin_plus_small, NULL, a, b, goal, tol, deg_limit, eval_limit, depth_limit, flags, prec);
                break;

            case 12:
                acb_set_d(a, -1020.0);
                acb_set_d(b, -1010.0);
                acb_calc_integrate(s, f_exp, NULL, a, b, goal, tol, deg_limit, eval_limit, depth_limit, flags, prec);
                break;

            case 13:
                acb_set_d(a, 0);
                acb_set_d(b, ceil(sqrt(goal * 0.693147181) + 1.0));
                acb_calc_integrate(s, f_gaussian, NULL, a, b, goal, tol, deg_limit, eval_limit, depth_limit, flags, prec);
                acb_mul(b, b, b, prec);
                acb_neg(b, b);
                acb_exp(b, b, prec);
                arb_add_error(acb_realref(s), acb_realref(b));
                break;

            case 14:
                acb_set_d(a, 0.0);
                acb_set_d(b, 1.0);
                acb_calc_integrate(s, f_wolfram, NULL, a, b, goal, tol, deg_limit, eval_limit, depth_limit, flags, prec);
                break;

            case 15:
                acb_set_d(a, 0.0);
                acb_set_d(b, 8.0);
                acb_calc_integrate(s, f_monster, NULL, a, b, goal, tol, deg_limit, eval_limit, depth_limit, flags, prec);
                break;

            case 16:
                acb_set_d(a, 0);
                acb_set_d(b, ceil(goal * 0.693147181 + 1.0));
                acb_calc_integrate(s, f_sech, NULL, a, b, goal, tol, deg_limit, eval_limit, depth_limit, flags, prec);
                acb_neg(b, b);
                acb_exp(b, b, prec);
                acb_mul_2exp_si(b, b, 1);
                arb_add_error(acb_realref(s), acb_realref(b));
                break;

            case 17:
                acb_set_d(a, 0);
                acb_set_d(b, ceil(goal * 0.693147181 / 3.0 + 2.0));
                acb_calc_integrate(s, f_sech3, NULL, a, b, goal, tol, deg_limit, eval_limit, depth_limit, flags, prec);
                acb_neg(b, b);
                acb_mul_ui(b, b, 3, prec);
                acb_exp(b, b, prec);
                acb_mul_2exp_si(b, b, 3);
                acb_div_ui(b, b, 3, prec);
                arb_add_error(acb_realref(s), acb_realref(b));
                break;

            default:
                abort();
            }
            TIMEIT_ONCE_STOP
        }
        flint_printf("I%d = ", integral);
        acb_printn(s, 3.333 * prec, 0);
        flint_printf("\n\n");
    }

    acb_clear(a);
    acb_clear(b);
    acb_clear(s);
    acb_clear(t);
    mag_clear(tol);

    flint_cleanup();
    return 0;
}

