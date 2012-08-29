#include "mpr.h"

int main()
{
    long iter;
    flint_rand_t state;

    printf("randtest....");
    fflush(stdout);

    flint_randinit(state);

    /* test exact roundtrip (todo: test rounding) */
    for (iter = 0; iter < 10000; iter++)
    {
        mpr_t x;

        mpr_init(x);
        mpr_randtest(x, state, 1 + n_randint(state, 10));

        if (!mpr_is_normalized(x))
        {
            printf("FAIL: not normalized\n");
            abort();
        }

        mpr_clear(x);
    }

    printf("PASS\n");
    flint_randclear(state);
    return EXIT_SUCCESS;
}
