#include "mpr.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("set_mpfr....");
    fflush(stdout);

    flint_randinit(state);
    _flint_rand_init_gmp(state);

    /* underscore methods: test exact roundtrip (todo: test rounding) */
    for (iter = 0; iter < 10000; iter++)
    {
        mp_ptr x;
        long exp;
        mpfr_t t, u;
        long n;
        mpfr_rnd_t rnd;

        n = 1 + n_randint(state, 10);
        x = malloc(n * sizeof(mp_limb_t));
        mpfr_init2(t, n * FLINT_BITS);
        mpfr_init2(u, n * FLINT_BITS);

        mpfr_urandom(t, state->gmp_state, MPFR_RNDN);
        mpfr_mul_2exp(t, t, n_randint(state, 20), MPFR_RNDN);
        mpfr_div_2exp(t, t, n_randint(state, 20), MPFR_RNDN);

        rnd = MPFR_RNDD;
        exp = _mpr_set_mpfr(x, t, n);
        _mpr_get_mpfr(u, x, exp, n, rnd);

        if (mpfr_cmp(t, u) != 0)
        {
            printf("FAIL!\n");
            mpfr_printf("\n%.200Rf\n", t);
            mpfr_printf("\n%.200Rf\n", u);
            abort();
        }

        mpfr_clear(t);
        mpfr_clear(u);
        free(x);
    }

    /* test exact roundtrip (todo: test rounding) */
    for (iter = 0; iter < 10000; iter++)
    {
        mpr_t x;
        mpfr_t t, u;
        long n;

        n = 1 + n_randint(state, 10);
        mpfr_init2(t, n * FLINT_BITS);
        mpfr_init2(u, n * FLINT_BITS);

        mpr_init(x);

        mpfr_urandom(t, state->gmp_state, MPFR_RNDN);
        mpfr_mul_2exp(t, t, n_randint(state, 20), MPFR_RNDN);
        mpfr_div_2exp(t, t, n_randint(state, 20), MPFR_RNDN);

        if (n_randint(state, 2))
            mpfr_neg(t, t, MPFR_RNDN);

        mpr_set_mpfr(x, t);

        if (!mpr_is_normalized(x))
        {
            printf("FAIL: not normalized\n");
            abort();
        }

        mpr_get_mpfr(u, x, MPFR_RNDN);

        if (mpfr_cmp(t, u) != 0)
        {
            printf("FAIL!\n");
            mpfr_printf("\n%.200Rf\n", t);
            mpfr_printf("\n%.200Rf\n", u);
            abort();
        }

        mpfr_clear(t);
        mpfr_clear(u);
        mpr_clear(x);
    }

    printf("PASS\n");

    flint_randclear(state);
    _fmpz_cleanup();
    return EXIT_SUCCESS;
}
