/* This file is public domain. Author: Fredrik Johansson. */

#include <string.h>
#include "fmprb_calc.h"
#include "profiler.h"

long eval_count = 0;

int
z_function(fmprb_ptr out, const fmprb_t inp, void * params, long order, long prec)
{
    fmprb_struct x[2];

    fmprb_init(x);
    fmprb_init(x + 1);

    fmprb_set(x, inp);
    fmprb_one(x + 1);

    _fmprb_poly_riemann_siegel_z_series(out, x, FLINT_MIN(2, order), order, prec);

    fmprb_clear(x);
    fmprb_clear(x + 1);

    eval_count++;
    return 0;
}

int
sin_x(fmprb_ptr out, const fmprb_t inp, void * params, long order, long prec)
{
    int xlen = FLINT_MIN(2, order);

    fmprb_set(out, inp);
    if (xlen > 1)
        fmprb_one(out + 1);

    _fmprb_poly_sin_series(out, out, xlen, order, prec);

    eval_count++;
    return 0;
}

int
sin_x2(fmprb_ptr out, const fmprb_t inp, void * params, long order, long prec)
{
    fmprb_ptr x;

    int xlen = FLINT_MIN(2, order);
    int ylen = FLINT_MIN(3, order);

    x = _fmprb_vec_init(xlen);

    fmprb_set(x, inp);
    if (xlen > 1)
        fmprb_one(x + 1);

    _fmprb_poly_mullow(out, x, xlen, x, xlen, ylen, prec);
    _fmprb_poly_sin_series(out, out, ylen, order, prec);

    _fmprb_vec_clear(x, xlen);

    eval_count++;
    return 0;
}

int
sin_1x(fmprb_ptr out, const fmprb_t inp, void * params, long order, long prec)
{
    fmprb_ptr x;
    int xlen = FLINT_MIN(2, order);

    x = _fmprb_vec_init(xlen);

    fmprb_set(x, inp);
    if (xlen > 1)
        fmprb_one(x + 1);

    _fmprb_poly_inv_series(out, x, xlen, order, prec);
    _fmprb_poly_sin_series(out, out, order, order, prec);

    _fmprb_vec_clear(x, xlen);

    eval_count++;
    return 0;
}

int main(int argc, char *argv[])
{
    fmprb_ptr blocks;
    fmprb_calc_func_t function;
    int * info;
    long digits, low_prec, high_prec, i, num, found_roots, found_unknown;
    long maxdepth, maxeval, maxfound;
    int refine;
    double a, b;
    fmprb_t t;
    fmpr_t fa, fb, C;

    if (argc < 4)
    {
        printf("real_roots function a b [-refine d] [-verbose] "
            "[-maxdepth n] [-maxeval n] [-maxfound n] [-prec n]\n");
        printf("available functions:\n");
        printf("  0  Z(x), Riemann-Siegel Z-function\n");
        printf("  1  sin(x)\n");
        printf("  2  sin(x^2)\n");
        printf("  3  sin(1/x)\n");
        return 1;
    }

    switch (atoi(argv[1]))
    {
        case 0:
            function = z_function;
            break;
        case 1:
            function = sin_x;
            break;
        case 2:
            function = sin_x2;
            break;
        case 3:
            function = sin_1x;
            break;
        default:
            printf("require a function 0-3\n");
            return 1;
    }

    a = atof(argv[2]);
    b = atof(argv[3]);

    if (a >= b)
    {
        printf("require a < b!\n");
        return 1;
    }

    refine = 0;
    digits = 0;
    maxdepth = 30;
    maxeval = 100000;
    maxfound = 100000;
    low_prec = 30;

    for (i = 4; i < argc; i++)
    {
        if (!strcmp(argv[i], "-refine"))
        {
            refine = 1;
            digits = atol(argv[i+1]);
        }
        else if (!strcmp(argv[i], "-verbose"))
        {
            fmprb_calc_verbose = 1;
        }
        else if (!strcmp(argv[i], "-maxdepth"))
        {
            maxdepth = atol(argv[i+1]);
        }
        else if (!strcmp(argv[i], "-maxeval"))
        {
            maxeval = atol(argv[i+1]);
        }
        else if (!strcmp(argv[i], "-maxfound"))
        {
            maxfound = atol(argv[i+1]);
        }
        else if (!strcmp(argv[i], "-prec"))
        {
            low_prec = atol(argv[i+1]);
        }
    }

    high_prec = digits * 3.32192809488736 + 10;
    found_roots = 0;
    found_unknown = 0;

    fmprb_init(t);
    fmpr_init(C);
    fmpr_init(fa);
    fmpr_init(fb);

    fmpr_set_d(fa, a);
    fmpr_set_d(fb, b);

    fmpr_add(fmprb_midref(t), fb, fa, FMPR_PREC_EXACT, FMPR_RND_DOWN);
    fmpr_sub(fmprb_radref(t), fb, fa, FMPR_PREC_EXACT, FMPR_RND_DOWN);
    fmprb_mul_2exp_si(t, t, -1);

    printf("interval: "); fmprb_printd(t, 15); printf("\n");
    printf("maxdepth = %ld, maxeval = %ld, maxfound = %ld, low_prec = %ld\n",
        maxdepth, maxeval, maxfound, low_prec);

    TIMEIT_ONCE_START

    num = fmprb_calc_isolate_roots(&blocks, &info, function,
        NULL, t, maxdepth, maxeval, maxfound, low_prec);

    for (i = 0; i < num; i++)
    {
        if (info[i] != 1)
        {
            if (fmprb_calc_verbose)
            {
                printf("unable to count roots in ");
                fmprb_printd(blocks + i, 15);
                printf("\n");
            }
            found_unknown++;
            continue;
        }

        found_roots++;

        if (!refine)
            continue;

        if (fmprb_calc_refine_root_bisect(t,
            function, NULL, blocks + i, 5, low_prec)
            != FMPRB_CALC_SUCCESS)
        {
            printf("warning: some bisection steps failed!\n");
        }

        if (fmprb_calc_verbose)
        {
            printf("after bisection 1: ");
            fmprb_printd(t, 15);
            printf("\n");
        }

        if (fmprb_calc_refine_root_bisect(blocks + i,
            function, NULL, t, 5, low_prec)
            != FMPRB_CALC_SUCCESS)
        {
            printf("warning: some bisection steps failed!\n");
        }

        if (fmprb_calc_verbose)
        {
            printf("after bisection 2: ");
            fmprb_printd(blocks + i, 15);
            printf("\n");
        }

        fmprb_calc_newton_conv_factor(C, function, NULL, t, low_prec);

        if (fmprb_calc_refine_root_newton(blocks + i, function, NULL,
            blocks + i, t, C, 10, high_prec) != FMPRB_CALC_SUCCESS)
        {
            printf("warning: some newton steps failed!\n");
        }

        printf("refined root:\n");
        fmprb_printd(blocks + i, digits + 2);
        printf("\n\n");
    }

    printf("---------------------------------------------------------------\n");
    printf("Found roots: %ld\n", found_roots);
    printf("Subintervals possibly containing undetected roots: %ld\n", found_unknown);
    printf("Function evaluations: %ld\n", eval_count);

    TIMEIT_ONCE_STOP
    SHOW_MEMORY_USAGE

    fmprb_clear(t);
    fmpr_clear(C);
    fmpr_clear(fa);
    fmpr_clear(fb);
    flint_cleanup();
    return 0;
}

