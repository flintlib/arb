#include "mpr.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("set_ui....");
    fflush(stdout);

    flint_randinit(state);
    _flint_rand_init_gmp(state);

    for (iter = 0; iter < 10000; iter++)
    {
        mpr_t x, y;
        mpfr_t t;
        ulong c;

        c = n_randtest(state);

        mpr_init(x);
        mpr_init(y);
        mpfr_init2(t, FLINT_BITS);

        mpr_set_ui(x, c);
        mpfr_set_ui(t, c, MPFR_RNDN);
        mpr_set_mpfr(y, t);

        if (!mpr_equal(x, y))
        {
            printf("FAIL!\n");
            printf("c = %lu\n", c);
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
