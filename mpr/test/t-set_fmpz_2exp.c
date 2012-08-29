#include "mpr.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("set_fmpz_2exp....");
    fflush(stdout);

    flint_randinit(state);
    _flint_rand_init_gmp(state);

    for (iter = 0; iter < 10000; iter++)
    {
        mpr_t x, y;
        mpfr_t t;
        fmpz_t c;
        mpz_t d;
        long exp;

        fmpz_init(c);
        mpz_init(d);

        fmpz_randtest(c, state, 200);
        fmpz_get_mpz(d, c);
        exp = n_randint(state, 100) - 50;

        mpr_init(x);
        mpr_init(y);
        mpfr_init2(t, 200);

        mpr_set_fmpz_2exp(x, c, exp);
        mpfr_set_z_2exp(t, d, exp, MPFR_RNDN);
        mpr_set_mpfr(y, t);

        if (!mpr_equal(x, y))
        {
            printf("FAIL!\n");
            printf("exp = %ld, c = ", exp); fmpz_print(c); printf("\n");
            printf("x = "); mpr_debug(x);
            printf("y = "); mpr_debug(y);
            abort();
        }

        mpr_clear(x);
        mpr_clear(y);
        mpfr_clear(t);

        fmpz_clear(c);
        mpz_clear(d);
    }

    printf("PASS\n");
    flint_randclear(state);
    _fmpz_cleanup();
    return EXIT_SUCCESS;
}
