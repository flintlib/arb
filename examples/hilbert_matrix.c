/* This file is public domain. Author: Fredrik Johansson. */

#include "fmprb_mat.h"
#include "profiler.h"

int main(int argc, char *argv[])
{
    fmprb_mat_t A;
    fmprb_t det;
    long i, j, prec, n;

    if (argc < 2)
    {
        printf("usage: build/examples/hilbert_matrix n\n");
        return 1;
    }

    n = atol(argv[1]);

    fmprb_mat_init(A, n, n);
    fmprb_init(det);

    TIMEIT_ONCE_START

    for (prec = 20; ; prec *= 2)
    {
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                fmprb_set_ui(fmprb_mat_entry(A, i, j), 1),
                fmprb_div_ui(fmprb_mat_entry(A, i, j),
                    fmprb_mat_entry(A, i, j), i + j + 1, prec);
            }
        }

        printf("prec=%ld: ", prec);

        fmprb_mat_det(det, A, prec);

        fmprb_printd(det, 10);
        printf("\n");

        if (!fmprb_contains_zero(det))
        {
            printf("success!\n");
            break;
        }
    }

    TIMEIT_ONCE_STOP

    SHOW_MEMORY_USAGE

    fmprb_mat_clear(A);
    fmprb_clear(det);
    flint_cleanup();
    return 0;
}

