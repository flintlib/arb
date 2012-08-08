#include "mprb.h"


/*
Let the input be [a-b, a+b] * 2^e. We require a > b >= 0 (otherwise the
interval contains zero or a negative number and the logarithm is not
defined). The error is largest at a-b, and we have

log(a * 2^e) - log((a-b) * 2^e) = log(1 + b/(a-b)).
*/

void
mprb_log_error(ufloat_t err, const mprb_t x)
{
    ufloat_t t, u;

    mprb_get_lower_bound_ufloat(t, x);
    ufloat_div(u, &x->rad, t);
    ufloat_log1p(err, u);
}

void
mprb_log_using_mpfr(mprb_t y, const mprb_t x)
{
    long prec;
    mpfr_t t, u;
    int input_approx, value_approx;
    ufloat_t err;

    if (mprb_contains_zero(x) || x->mid.sign == MPRB_SIGN_MINUS)
    {
        printf("mprb_log: interval contains zero or negative numbers\n");
        abort();
    }

    prec = y->mid.alloc * FLINT_BITS;

    mpfr_init2(t, x->mid.size * FLINT_BITS);
    mpfr_init2(u, prec);

    mprb_get_mid_mpfr(t, x, MPFR_RNDN);  /* exact */

    input_approx = !mprb_is_exact(x);

    if (input_approx)
        mprb_log_error(err, x);

    /* todo: handle exact case (log(1)) */
    mpfr_log(u, t, MPFR_RNDN);

    mprb_set_mpfr(y, u);

    ufloat_set_2exp(&y->rad, mprb_ulp_exp(y));

    if (input_approx)
        ufloat_add(&y->rad, &y->rad, err);

    mpfr_clear(t);
    mpfr_clear(u);
}
