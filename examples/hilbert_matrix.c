/* This file is public domain. Author: Fredrik Johansson. */

#include "arb_mat.h"
#include "profiler.h"

int main(int argc, char *argv[])
{
    arb_mat_t A;
    arb_t det;
    long i, j, prec, n;

    if (argc < 2)
    {
        printf("usage: build/examples/hilbert_matrix n\n");
        return 1;
    }

    n = atol(argv[1]);

    arb_mat_init(A, n, n);
    arb_init(det);

    TIMEIT_ONCE_START

    for (prec = 20; ; prec *= 2)
    {
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                arb_set_ui(arb_mat_entry(A, i, j), 1),
                arb_div_ui(arb_mat_entry(A, i, j),
                    arb_mat_entry(A, i, j), i + j + 1, prec);
            }
        }

        printf("prec=%ld: ", prec);

        arb_mat_det(det, A, prec);

        arb_printd(det, 10);
        printf("\n");

        if (!arb_contains_zero(det))
        {
            printf("success!\n");
            break;
        }
    }

    TIMEIT_ONCE_STOP

    SHOW_MEMORY_USAGE

    arb_mat_clear(A);
    arb_clear(det);
    flint_cleanup();
    return 0;
}

