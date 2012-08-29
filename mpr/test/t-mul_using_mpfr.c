#include "mpr.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("mul_using_mpfr....");
    fflush(stdout);

    flint_randinit(state);
    _flint_rand_init_gmp(state);

    for (iter = 0; iter < 100000; iter++)
    {
        mpr_t x, y, z;
        mpfr_t t, u, v, w, v_exact, error;
        long b1, b2, b3, error_exp;
        mpfr_rnd_t rnd;

        b1 = 2 + n_randint(state, 200);
        b2 = 2 + n_randint(state, 200);
        b3 = 2 + n_randint(state, 200);

        mpfr_init2(t, b1);
        mpfr_init2(u, b2);
        mpfr_init2(v, b3);
        mpfr_init2(w, b3);
        mpfr_init2(v_exact, b3);
        mpfr_init2(error, b3);

        mpr_init(x);
        mpr_init(y);
        mpr_init(z);

        mpfr_urandom(t, state->gmp_state, MPFR_RNDN);
        mpfr_mul_2exp(t, t, n_randint(state, 20), MPFR_RNDN);
        mpfr_div_2exp(t, t, n_randint(state, 20), MPFR_RNDN);
        if (n_randint(state, 2))
            mpfr_neg(u, u, MPFR_RNDN);

        mpfr_urandom(u, state->gmp_state, MPFR_RNDN);
        mpfr_mul_2exp(u, u, n_randint(state, 20), MPFR_RNDN);
        mpfr_div_2exp(u, u, n_randint(state, 20), MPFR_RNDN);
        if (n_randint(state, 2))
            mpfr_neg(u, u, MPFR_RNDN);

        mpr_randtest(z, state, 1 + n_randint(state, 200));

        switch (n_randint(state, 5))
        {
            case 0: rnd = MPFR_RNDN; break;
            case 1: rnd = MPFR_RNDD; break;
            case 2: rnd = MPFR_RNDU; break;
            case 3: rnd = MPFR_RNDA; break;
            default: rnd = MPFR_RNDZ;
        }

        mpr_set_mpfr(x, t);
        mpr_set_mpfr(y, u);

        error_exp = mpr_mul_using_mpfr(z, x, y, b3, rnd);

        if (mpfr_mul(v, t, u, rnd) == 0)
        {
            if (error_exp != LONG_MIN)
            {
                printf("FAIL: exact mpfr addition but error_exp != LONG_MIN\n");
                abort();
            }
        }
        else
        {
            /* check error bound */
            while (mpfr_mul(v_exact, t, u, rnd) != 0)
                mpfr_set_prec(v_exact, mpfr_get_prec(v_exact) * 2);

            while (mpfr_sub(error, v_exact, v, rnd) != 0)
                mpfr_set_prec(error, mpfr_get_prec(error) * 2);

            if (mpfr_get_exp(error) > error_exp)
            {
                printf("FAIL: exp(error) = %ld, error_exp = %ld\n",
                    mpfr_get_exp(error), error_exp);
                abort();
            }
        }

        mpr_get_mpfr(w, z, MPFR_RNDN);

        if (!mpr_is_normalized(z))
        {
            printf("FAIL: not normalized\n");
            abort();
        }

        if (mpfr_cmp(v, w) != 0)
        {
            printf("FAIL!\n");
            mpfr_printf("\n%.200Rf\n", t);
            mpfr_printf("\n%.200Rf\n", u);
            mpfr_printf("\n%.200Rf\n", v);
            mpfr_printf("\n%.200Rf\n", w);
            abort();
        }

        mpfr_clear(t);
        mpfr_clear(u);
        mpfr_clear(v);
        mpfr_clear(w);
        mpfr_clear(v_exact);
        mpfr_clear(error);
        mpr_clear(x);
        mpr_clear(y);
        mpr_clear(z);
    }

    printf("PASS\n");
    flint_randclear(state);
    return EXIT_SUCCESS;
}
