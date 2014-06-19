/* This file is public domain. Author: Fredrik Johansson. */

#include "fmpcb_calc.h"
#include "profiler.h"

int
sinx(fmpcb_ptr out, const fmpcb_t inp, void * params, long order, long prec)
{
    int xlen = FLINT_MIN(2, order);
    fmpcb_set(out, inp);
    if (xlen > 1)
        fmpcb_one(out + 1);
    _fmpcb_poly_sin_series(out, out, xlen, order, prec);
    return 0;
}

int
elliptic(fmpcb_ptr out, const fmpcb_t inp, void * params, long order, long prec)
{
    fmpcb_ptr t;
    t = _fmpcb_vec_init(order);
    fmpcb_set(t, inp);
    if (order > 1)
        fmpcb_one(t + 1);
    _fmpcb_poly_sin_series(t, t, FLINT_MIN(2, order), order, prec);
    _fmpcb_poly_mullow(out, t, order, t, order, order, prec);
    _fmpcb_vec_scalar_mul_2exp_si(t, out, order, -1);
    fmpcb_sub_ui(t, t, 1, prec);
    _fmpcb_vec_neg(t, t, order);
    _fmpcb_poly_rsqrt_series(out, t, order, order, prec);
    _fmpcb_vec_clear(t, order);
    return 0;
}

int
bessel(fmpcb_ptr out, const fmpcb_t inp, void * params, long order, long prec)
{
    fmpcb_ptr t;
    fmpcb_t z;
    ulong n;

    t = _fmpcb_vec_init(order);
    fmpcb_init(z);

    fmpcb_set(t, inp);
    if (order > 1)
        fmpcb_one(t + 1);

    n = 10;
    fmprb_set_si(fmpcb_realref(z), 20);
    fmprb_set_si(fmpcb_imagref(z), 10);

    /* z sin(t) */
    _fmpcb_poly_sin_series(out, t, FLINT_MIN(2, order), order, prec);
    _fmpcb_vec_scalar_mul(out, out, order, z, prec);

    /* t n */
    _fmpcb_vec_scalar_mul_ui(t, t, FLINT_MIN(2, order), n, prec);

    _fmpcb_poly_sub(out, t, FLINT_MIN(2, order), out, order, prec);

    _fmpcb_poly_cos_series(out, out, order, order, prec);

    _fmpcb_vec_clear(t, order);
    fmpcb_clear(z);
    return 0;
}

int main(int argc, char *argv[])
{
    fmpcb_t r, s, a, b;
    fmpr_t inr, outr;
    long digits, prec;

    if (argc < 2)
    {
        printf("integrals d\n");
        printf("compute integrals using d decimal digits of precision\n");
        return 1;
    }

    fmpcb_init(r);
    fmpcb_init(s);
    fmpcb_init(a);
    fmpcb_init(b);
    fmpr_init(inr);
    fmpr_init(outr);

    fmprb_calc_verbose = 0;

    digits = atol(argv[1]);
    prec = digits * 3.32193;
    printf("Digits: %ld\n", digits);

    printf("----------------------------------------------------------------\n");
    printf("Integral of sin(t) from 0 to 100.\n");
    fmpr_set_d(inr, 0.125);
    fmpr_set_d(outr, 1.0);
    TIMEIT_ONCE_START
    fmpcb_set_si(a, 0);
    fmpcb_set_si(b, 100);
    fmpcb_calc_integrate_taylor(r, sinx, NULL, a, b, inr, outr, prec, 1.1 * prec);
    printf("RESULT:\n");
    fmpcb_printd(r, digits); printf("\n");
    TIMEIT_ONCE_STOP

    printf("----------------------------------------------------------------\n");
    printf("Elliptic integral F(phi, m) = integral of 1/sqrt(1 - m*sin(t)^2)\n");
    printf("from 0 to phi, with phi = 6+6i, m = 1/2. Integration path\n");
    printf("0 -> 6 -> 6+6i.\n");
    fmpr_set_d(inr, 0.2);
    fmpr_set_d(outr, 0.5);
    TIMEIT_ONCE_START
    fmpcb_set_si(a, 0);
    fmpcb_set_si(b, 6);
    fmpcb_calc_integrate_taylor(r, elliptic, NULL, a, b, inr, outr, prec, 1.1 * prec);
    fmpcb_set_si(a, 6);
    fmprb_set_si(fmpcb_realref(b), 6);
    fmprb_set_si(fmpcb_imagref(b), 6);
    fmpcb_calc_integrate_taylor(s, elliptic, NULL, a, b, inr, outr, prec, 1.1 * prec);
    fmpcb_add(r, r, s, prec);
    printf("RESULT:\n");
    fmpcb_printd(r, digits); printf("\n");
    TIMEIT_ONCE_STOP

    printf("----------------------------------------------------------------\n");
    printf("Bessel function J_n(z) = (1/pi) * integral of cos(t*n - z*sin(t))\n");
    printf("from 0 to pi. With n = 10, z = 20 + 10i.\n");
    fmpr_set_d(inr, 0.1);
    fmpr_set_d(outr, 0.5);
    TIMEIT_ONCE_START
    fmpcb_set_si(a, 0);
    fmpcb_const_pi(b, 3 * prec);
    fmpcb_calc_integrate_taylor(r, bessel, NULL, a, b, inr, outr, prec, 3 * prec);
    fmpcb_div(r, r, b, prec);
    printf("RESULT:\n");
    fmpcb_printd(r, digits); printf("\n");
    TIMEIT_ONCE_STOP

    fmpcb_clear(r);
    fmpcb_clear(s);
    fmpcb_clear(a);
    fmpcb_clear(b);
    fmpr_clear(inr);
    fmpr_clear(outr);

    flint_cleanup();
    return 0;
}

