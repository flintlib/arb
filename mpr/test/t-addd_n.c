#include "mpr.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("addd_n....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000000; iter++)
    {
        mp_ptr x, y, z;
        mpfr_t t, u, v, w;
        long n, xexp, yexp, zexp;

        n = 1 + n_randint(state, 40);

        x = malloc(n * sizeof(mp_limb_t));
        y = malloc(n * sizeof(mp_limb_t));
        z = malloc(n * sizeof(mp_limb_t));

        mpfr_init2(t, n * FLINT_BITS);
        mpfr_init2(u, n * FLINT_BITS);
        mpfr_init2(v, n * FLINT_BITS);
        mpfr_init2(w, n * FLINT_BITS);

        _mpr_randtest(x, state, n);
        _mpr_randtest(y, state, n);
        _mpr_randtest(z, state, n);

        xexp = n_randint(state, 10) - 20;
        yexp = n_randint(state, 10) - 20;

        _mpr_get_mpfr(t, x, xexp, n, MPFR_RNDN);
        _mpr_get_mpfr(u, y, yexp, n, MPFR_RNDN);

        if (xexp >= yexp)
            zexp = xexp + _mpr_addd_n(z, x, y, xexp - yexp, n);
        else
            zexp = yexp + _mpr_addd_n(z, y, x, yexp - xexp, n);

        _mpr_get_mpfr(w, z, zexp, n, MPFR_RNDN);

        mpfr_add(v, t, u, MPFR_RNDD);

        if (mpfr_cmp(v, w) != 0)
        {
            printf("FAIL!\n");
            _mpr_printr(x, n);
            _mpr_printr(y, n);
            _mpr_printr(z, n);

            mpfr_printf("\n%.80Rf\n", t);
            mpfr_printf("\n%.80Rf\n", u);
            mpfr_printf("\n%.80Rf\n", v);
            mpfr_printf("\n%.80Rf\n", w);
            abort();
        }

        mpfr_clear(t);
        mpfr_clear(u);
        mpfr_clear(v);
        mpfr_clear(w);
        free(x);
        free(y);
        free(z);
    }

    printf("PASS\n");

    flint_randclear(state);
    return EXIT_SUCCESS;
}
