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

    printf("mul....");
    fflush(stdout);

    flint_randinit(state);
    _flint_rand_init_gmp(state);

    for (iter = 0; iter < 1000000; iter++)
    {
        mprb_t x, y, z;
        mpfr_t t, u, v;

        mprb_init(x, 1 + n_randint(state, 300));
        mprb_init(y, 1 + n_randint(state, 300));
        mprb_init(z, 1 + n_randint(state, 300));

        do { mprb_randtest(x, state, -10, 10); } while (x->size == 0);
        do { mprb_randtest(y, state, -10, 10); } while (y->size == 0);
        do { mprb_randtest(z, state, -10, 10); } while (z->size == 0);

        mpfr_init2(t, 1000);
        mpfr_init2(u, 1000);
        mpfr_init2(v, 1000);

        if (n_randint(state, 2))
            mprb_get_interval_mpfr(t, v, x);
        else
            mprb_get_interval_mpfr(v, t, x);

        if (n_randint(state, 2))
            mprb_get_interval_mpfr(u, v, y);
        else
            mprb_get_interval_mpfr(v, u, y);

        mprb_mul(z, x, y);
        mpfr_mul(v, t, u, MPFR_RNDN);

        if (!mprb_contains_mpfr(z, v))
        {
            printf("FAIL!\n");
            mprb_debug(x);
            mprb_debug(y);
            mprb_debug(z);

            mpfr_debug(t);
            mpfr_debug(u);
            mpfr_debug(v);

            abort();
        }

        mpfr_clear(t);
        mpfr_clear(u);
        mpfr_clear(v);

        mprb_clear(x);
        mprb_clear(y);
        mprb_clear(z);
    }

    printf("PASS\n");
    flint_randclear(state);
    _fmpz_cleanup();
    return EXIT_SUCCESS;
}