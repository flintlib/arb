/*
    Copyright (C) 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "partitions.h"
#include "arb.h"
#include "math.h"

#define DOUBLE_CUTOFF 40
#define DOUBLE_ERR 1e-12

#define DOUBLE_PREC 53
#define MIN_PREC 20
#define PI 3.141592653589793238462643
#define INV_LOG2 (1.44269504088896340735992468 + 1e-12)
#define HRR_A (1.1143183348516376904 + 1e-12)  /* 44*pi^2/(225*sqrt(3)) */
#define HRR_B (0.0592384391754448833 + 1e-12)  /* pi*sqrt(2)/75 */
#define HRR_C (2.5650996603237281911 + 1e-12)  /* pi*sqrt(2/3) */
#define HRR_D (1.2424533248940001551 + 1e-12)  /* log(2) + log(3)/2 */

static double
partitions_remainder_bound(double n, double terms)
{
    return HRR_A/sqrt(terms)
            + HRR_B*sqrt(terms/(n-1)) * sinh(HRR_C * sqrt(n)/terms);
}

/* Crude upper bound, sufficient to estimate the precision */
static double
log_sinh(double x)
{
    if (x > 4)
        return x;
    else
        return log(x) + x*x*(1/6.);
}

static double
partitions_remainder_bound_log2(double n, double N)
{
    double t1, t2;

    t1 = log(HRR_A) - 0.5*log(N);
    t2 = log(HRR_B) + 0.5*(log(N) - log(n-1)) + log_sinh(HRR_C * sqrt(n)/N);

    return (FLINT_MAX(t1, t2) + 1) * INV_LOG2;
}

slong
partitions_hrr_needed_terms(double n)
{
    slong N;
    for (N = 1; partitions_remainder_bound_log2(n, N) > 10; N++);
    for ( ; partitions_remainder_bound(n, N) > 0.4; N++);
    return N;
}

static double
partitions_term_bound(double n, double k)
{
    return ((PI*sqrt(24*n-1) / (6.0*k)) + HRR_D - log(24.0*n-1) + 0.5*log(k)) * INV_LOG2;
}

/* Bound number of prime factors in k */
static mp_limb_t primorial_tab[] = {
    1, 2, 6, 30, 210, 2310, 30030, 510510, 9699690, 223092870,
#if FLINT64
    UWORD(6469693230), UWORD(200560490130), UWORD(7420738134810), 
    UWORD(304250263527210), UWORD(13082761331670030), UWORD(614889782588491410)
#endif
};

static __inline__ int
bound_primes(ulong k)
{
    int i;

    for (i = 0; i < sizeof(primorial_tab) / sizeof(mp_limb_t); i++)
        if (k <= primorial_tab[i])
            return i;

    return i;
}


static __inline__ slong
log2_ceil(double x)
{
    /* ceil(log2(n)) = bitcount(n-1);
       this is too large if x is a power of two */
    return FLINT_BIT_COUNT((slong) x);
}

static slong
partitions_prec_bound(double n, slong k, slong N)
{
    slong prec;

    prec = partitions_term_bound(n, k);
    prec += log2_ceil(8 * N * (26 * (sqrt(n) / k) + 7 * bound_primes(k) + 22));

    return prec;
}

static double
cos_pi_pq(mp_limb_signed_t p, mp_limb_signed_t q)
{
    /* Force 0 <= p < q */
    p = FLINT_ABS(p);
    p %= (2 * q);
    if (p >= q)
        p = 2 * q - p;

    if (4 * p <= q)
        return cos(p * PI / q);
    else if (4 * p < 3 * q)
        return sin((q - 2*p) * PI / (2 * q));
    else
        return -cos((q - p) * PI / q);
}

static double
eval_trig_prod_d(trig_prod_t prod)
{
    int i;
    double s;
    mp_limb_t v;

    if (prod->prefactor == 0)
        return 0.0;

    s = prod->prefactor;

    v = n_gcd(FLINT_MAX(prod->sqrt_p, prod->sqrt_q),
              FLINT_MIN(prod->sqrt_p, prod->sqrt_q));
    prod->sqrt_p /= v;
    prod->sqrt_q /= v;

    if (prod->sqrt_p != 1)
        s *= sqrt(prod->sqrt_p);
    if (prod->sqrt_q != 1)
        s /= sqrt(prod->sqrt_q);

    for (i = 0; i < prod->n; i++)
        s *= cos_pi_pq(prod->cos_p[i], prod->cos_q[i]);

    return s;
}

static void
eval_trig_prod(arb_t sum, trig_prod_t prod, slong prec)
{
    int i;
    mp_limb_t v;
    arb_t t;

    if (prod->prefactor == 0)
    {
        arb_zero(sum);
        return;
    }

    arb_init(t);

    arb_set_si(sum, prod->prefactor);
    v = n_gcd(FLINT_MAX(prod->sqrt_p, prod->sqrt_q),
              FLINT_MIN(prod->sqrt_p, prod->sqrt_q));
    prod->sqrt_p /= v;
    prod->sqrt_q /= v;

    if (prod->sqrt_p != 1)
    {
        arb_sqrt_ui(t, prod->sqrt_p, prec);
        arb_mul(sum, sum, t, prec);
    }

    if (prod->sqrt_q != 1)
    {
        arb_rsqrt_ui(t, prod->sqrt_q, prec);
        arb_mul(sum, sum, t, prec);
    }

    for (i = 0; i < prod->n; i++)
    {
        fmpq_t pq;
        *fmpq_numref(pq) = prod->cos_p[i];
        *fmpq_denref(pq) = prod->cos_q[i];
        arb_cos_pi_fmpq(t, pq, prec);
        arb_mul(sum, sum, t, prec);
    }

    arb_clear(t);
}

static void
sinh_cosh_divk_precomp(arb_t sh, arb_t ch, arb_t ex, slong k, slong prec)
{
    arb_t t;
    arb_init(t);
    arb_set_round(t, ex, prec);
    arb_root_ui(ch, t, k, prec);
    /* The second term doesn't need full precision,
       but this doesn't affect performance that much... */
    arb_inv(t, ch, prec);
    arb_sub(sh, ch, t, prec);
    arb_add(ch, ch, t, prec);
    arb_mul_2exp_si(ch, ch, -1);
    arb_mul_2exp_si(sh, sh, -1);
    arb_clear(t);
}


void
partitions_hrr_sum_arb(arb_t x, const fmpz_t n, slong N0, slong N, int use_doubles)
{
    trig_prod_t prod;
    arb_t acc, C, t1, t2, t3, t4, exp1;
    fmpz_t n24;
    slong k, prec, res_prec, acc_prec, guard_bits;
    double nd, Cd;

    if (fmpz_cmp_ui(n, 2) <= 0)
    {
        flint_abort();
    }

    nd = fmpz_get_d(n);

    /* compute initial precision */
    guard_bits = 2 * FLINT_BIT_COUNT(N) + 32;
    prec = partitions_remainder_bound_log2(nd, N0) + guard_bits;
    prec = FLINT_MAX(prec, DOUBLE_PREC);
    res_prec = acc_prec = prec;

    arb_init(acc);
    arb_init(C);
    arb_init(t1);
    arb_init(t2);
    arb_init(t3);
    arb_init(t4);
    arb_init(exp1);
    fmpz_init(n24);

    arb_zero(x);

    /* n24 = 24n - 1 */
    fmpz_set(n24, n);
    fmpz_mul_ui(n24, n24, 24);
    fmpz_sub_ui(n24, n24, 1);

    /* C = (pi/6) sqrt(24n-1) */
    arb_const_pi(t1, prec);
    arb_sqrt_fmpz(t2, n24, prec);
    arb_mul(t1, t1, t2, prec);
    arb_div_ui(C, t1, 6, prec);

    /* exp1 = exp(C) */
    arb_exp(exp1, C, prec);

    Cd = PI * sqrt(24*nd-1) / 6;

    for (k = N0; k <= N; k++)
    {
        trig_prod_init(prod);
        arith_hrr_expsum_factored(prod, k, fmpz_fdiv_ui(n, k));

        if (prod->prefactor != 0)
        {

            if (prec > MIN_PREC)
                prec = partitions_prec_bound(nd, k, N);

            prod->prefactor *= 4;
            prod->sqrt_p *= 3;
            prod->sqrt_q *= k;

            if (prec > DOUBLE_CUTOFF || !use_doubles)
            {
                /* Compute A_k(n) * sqrt(3/k) * 4 / (24*n-1) */
                eval_trig_prod(t1, prod, prec);
                arb_div_fmpz(t1, t1, n24, prec);

                /* Multiply by (cosh(z) - sinh(z)/z) where z = C / k */
                arb_set_round(t2, C, prec);
                arb_div_ui(t2, t2, k, prec);

                if (k < 35 && prec > 1000)
                    sinh_cosh_divk_precomp(t3, t4, exp1, k, prec);
                else
                    arb_sinh_cosh(t3, t4, t2, prec);

                arb_div(t3, t3, t2, prec);
                arb_sub(t2, t4, t3, prec);
                arb_mul(t1, t1, t2, prec);
            }
            else
            {
                double xx, zz, xxerr;

                xx = eval_trig_prod_d(prod) / (24*nd - 1);
                zz = Cd / k;
                xx = xx * (cosh(zz) - sinh(zz) / zz);

                xxerr = fabs(xx) * DOUBLE_ERR + DOUBLE_ERR;
                arf_set_d(arb_midref(t1), xx);
                mag_set_d(arb_radref(t1), xxerr);
            }

            /* Add to accumulator */
            arb_add(acc, acc, t1, acc_prec);

            if (acc_prec > 2 * prec + 32)
            {
                arb_add(x, x, acc, res_prec);
                arb_zero(acc);
                acc_prec = prec + 32;
            }
        }
    }

    arb_add(x, x, acc, res_prec);

    fmpz_clear(n24);
    arb_clear(acc);
    arb_clear(exp1);
    arb_clear(C);
    arb_clear(t1);
    arb_clear(t2);
    arb_clear(t3);
    arb_clear(t4);
}


