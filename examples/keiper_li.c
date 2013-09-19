#include <string.h>
#include "fmprb.h"
#include "fmpcb.h"
#include "fmprb_poly.h"
#include "zeta.h"
#include "profiler.h"

void
keiper_li_series(fmprb_ptr z, long len, long prec)
{
    fmprb_ptr t, u, v;

    t = _fmprb_vec_init(len);
    u = _fmprb_vec_init(len);
    v = _fmprb_vec_init(len);

    /* -zeta(s) */
    printf("zeta: ");
    TIMEIT_ONCE_START
    fmprb_zero(t + 0);
    fmprb_one(t + 1);
    fmprb_one(u);
    _fmprb_poly_zeta_series(v, t, 2, u, 0, len, prec);
    _fmprb_vec_neg(v, v, len);
    TIMEIT_ONCE_STOP

    SHOW_MEMORY_USAGE

    /* logarithm */
    printf("log: ");
    TIMEIT_ONCE_START
    _fmprb_poly_log_series(t, v, len, len, prec);
    TIMEIT_ONCE_STOP

    /* add log(gamma(1+s/2)) */
    printf("gamma: ");
    TIMEIT_ONCE_START
    fmprb_one(u);
    fmprb_one(u + 1);
    fmprb_mul_2exp_si(u + 1, u + 1, -1);
    _fmprb_poly_lgamma_series(v, u, 2, len, prec);
    _fmprb_vec_add(t, t, v, len, prec);
    TIMEIT_ONCE_STOP

    /* subtract 0.5 s log(pi) */
    fmprb_const_pi(u, prec);
    fmprb_log(u, u, prec);
    fmprb_mul_2exp_si(u, u, -1);
    fmprb_sub(t + 1, t + 1, u, prec);

    /* add log(1-s) */
    fmprb_one(u);
    fmprb_set_si(u + 1, -1);
    _fmprb_poly_log_series(v, u, 2, len, prec);
    _fmprb_vec_add(t, t, v, len, prec);

    /* binomial transform */
    printf("binomial transform: ");
    TIMEIT_ONCE_START
    fmprb_set(z, t);
    _fmprb_vec_neg(t + 1, t + 1, len - 1);
    _fmprb_poly_binomial_transform(z + 1, t + 1, len - 1, len - 1, prec);
    TIMEIT_ONCE_STOP

    _fmprb_vec_clear(t, len);
    _fmprb_vec_clear(u, len);
    _fmprb_vec_clear(v, len);
}

int main(int argc, char *argv[])
{
    long i, len, prec, num_threads;
    char * out_file;
    fmprb_ptr z;

    if (argc < 2)
    {
        printf("keiper_li n [-prec prec] [-threads num_threads] [-out out_file]\n");
        return 1;
    }

    len = atol(argv[1]) + 1;
    prec = 1.1 * len + 50;
    num_threads = 1;
    out_file = NULL;

    for (i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-prec"))
            prec = atol(argv[i+1]);
        else if (!strcmp(argv[i], "-threads"))
            num_threads = atol(argv[i+1]);
        else if (!strcmp(argv[i], "-out"))
            out_file = argv[i+1];
    }

    flint_set_num_threads(num_threads);

    z = _fmprb_vec_init(len);

    keiper_li_series(z, len, prec);

    for (i = 0; i < len; i++)
    {
        if (i <= 10 || len - i <= 10)
        {
            printf("%ld: ", i); fmprb_printd(z + i, 50); printf("\n");
        }
    }

    SHOW_MEMORY_USAGE

    if (out_file != NULL)
    {
        FILE * fp = fopen(out_file, "w");
        for (i = 0; i < len; i++)
        {
            fprintf(fp, "%ld ", i);
            fmpz_fprint(fp, fmpr_manref(fmprb_midref(z + i)));
            fprintf(fp, " ");
            fmpz_fprint(fp, fmpr_expref(fmprb_midref(z + i)));
            fprintf(fp, " ");
            fmpz_fprint(fp, fmpr_manref(fmprb_radref(z + i)));
            fprintf(fp, " ");
            fmpz_fprint(fp, fmpr_expref(fmprb_radref(z + i)));
            fprintf(fp, "\n");
        }
        fclose(fp);
    }

    _fmprb_vec_clear(z, len);
    flint_cleanup();
    return 0;
}

