#include "mpr.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("mul_exact....");
    fflush(stdout);

    flint_randinit(state);
    _flint_rand_init_gmp(state);

    for (iter = 0; iter < 100000; iter++)
    {
        mpr_t x, y, z;
        mpfr_t t, u, v, w;
        long m, n;

        m = 1 + n_randint(state, 30);
        n = 1 + n_randint(state, 30);

        mpfr_init2(t, m * FLINT_BITS);
        mpfr_init2(u, n * FLINT_BITS);
        mpfr_init2(v, (m + n) * FLINT_BITS);
        mpfr_init2(w, (m + n) * FLINT_BITS);

        mpr_init(x);
        mpr_init(y);
        mpr_init(z);

        mpr_randtest(x, state, m);
        mpr_randtest(y, state, n);
        mpr_randtest(z, state, 1 + n_randint(state, 30));

        mpr_mul_exact(z, x, y);

        mpr_get_mpfr(t, x, MPFR_RNDN);
        mpr_get_mpfr(u, y, MPFR_RNDN);

        mpfr_mul(v, t, u, MPFR_RNDN);

        mpr_get_mpfr(w, z, MPFR_RNDN);

        if (!mpr_is_normalized(z))
        {
            printf("FAIL: not normalized\n");
            abort();
        }

        if (mpfr_cmp(v, w) != 0)
        {
            printf("FAIL!\n");
            mpfr_printf("\n%.500Rf\n", v);
            mpfr_printf("\n%.500Rf\n", w);
            abort();
        }

        mpfr_clear(t);
        mpfr_clear(u);
        mpfr_clear(v);
        mpfr_clear(w);
        mpr_clear(x);
        mpr_clear(y);
        mpr_clear(z);
    }

    printf("PASS\n");

    flint_randclear(state);
    _fmpz_cleanup();
    return EXIT_SUCCESS;
}
