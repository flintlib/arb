#include "mprb.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("set_mpfr....");
    fflush(stdout);

    flint_randinit(state);
    _flint_rand_init_gmp(state);

    for (iter = 0; iter < 10000; iter++)
    {
        mprb_t x;
        mpfr_t t;

        mprb_init(x, 1 + n_randint(state, 400));
        /* mprb_randtest(x, state); */

        mpfr_init2(t, 2 + n_randint(state, 400));

        mpfr_urandom(t, state->gmp_state, MPFR_RNDN);
        if (n_randint(state, 2))
            mpfr_neg(t, t, MPFR_RNDN);

        mpfr_mul_2exp(t, t, n_randint(state, 20), MPFR_RNDN);
        mpfr_div_2exp(t, t, n_randint(state, 20), MPFR_RNDN);

        mprb_set_mpfr(x, t);

        if (!mprb_contains_mpfr(x, t))
        {
            printf("FAIL!\n");
            mprb_debug(x);
            mpfr_printf("\n%Rf\n", t);
            mpn_debug(t->_mpfr_d, (t->_mpfr_prec + FLINT_BITS - 1) / FLINT_BITS);
            abort();
        }

        mpfr_clear(t);
        mprb_clear(x);
    }

    printf("PASS\n");
    flint_randclear(state);
    _fmpz_cleanup();
    return EXIT_SUCCESS;
}
