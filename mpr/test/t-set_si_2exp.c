#include "mpr.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("set_si_2exp....");
    fflush(stdout);

    flint_randinit(state);
    _flint_rand_init_gmp(state);

    for (iter = 0; iter < 10000; iter++)
    {
        mpr_t x, y;
        mpfr_t t;
        long c, exp;

        c = n_randtest(state);
        exp = n_randint(state, 100) - 50;

        mpr_init(x);
        mpr_init(y);
        mpfr_init2(t, FLINT_BITS);

        mpr_set_si_2exp(x, c, exp);
        mpfr_set_si_2exp(t, c, exp, MPFR_RNDN);
        mpr_set_mpfr(y, t);

        if (!mpr_equal(x, y))
        {
            printf("FAIL!\n");
            printf("c = %ld, exp = %ld\n", c, exp);
            printf("x = "); mpr_debug(x);
            printf("y = "); mpr_debug(y);
            abort();
        }

        mpr_clear(x);
        mpr_clear(y);
        mpfr_clear(t);
    }

    printf("PASS\n");
    flint_randclear(state);
    _fmpz_cleanup();
    return EXIT_SUCCESS;
}
