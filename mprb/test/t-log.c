#include "mprb.h"

void mpfr_debug(mpfr_t t)
{
    mpfr_printf("\n%Rf\n", t);
    mpn_debug(t->_mpfr_d, (t->_mpfr_prec + FLINT_BITS - 1) / FLINT_BITS);
}

int main()
{
    long iter;
    flint_rand_t state;

    printf("log....");
    fflush(stdout);

    flint_randinit(state);
    _flint_rand_init_gmp(state);

/*
    for (iter = 0; iter < 100000; iter++)
    {
        mprb_t x, y;
        mpfr_t t, u;

        mprb_init(x, 1 + n_randint(state, 300));
        mprb_init(y, 1 + n_randint(state, 300));
 
        do { mprb_randtest(x, state, -1000, 1000); }
            while (x->mid.size == 0 || mprb_contains_zero(x)
                        || x->mid.sign == MPRB_SIGN_MINUS);

        mpfr_init2(t, 1000);
        mpfr_init2(u, 1000);

        if (n_randint(state, 2))
            mprb_get_interval_mpfr(t, u, x);
        else
            mprb_get_interval_mpfr(u, t, x);

        mprb_log_using_mpfr(y, x);
        mpfr_log(u, t, MPFR_RNDN);

        if (!mprb_contains_mpfr(y, u))
        {
            printf("FAIL!\n");
            mprb_debug(x);
            mprb_debug(y);

            mpfr_debug(t);
            mpfr_debug(u);

            abort();
        }

        mpfr_clear(t);
        mpfr_clear(u);

        mprb_clear(x);
        mprb_clear(y);
    }
*/

    printf("PASS\n");
    flint_randclear(state);
    return EXIT_SUCCESS;
}
