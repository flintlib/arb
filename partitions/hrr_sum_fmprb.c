/*=============================================================================

    This file is part of ARB.

    ARB is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    ARB is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ARB; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Fredrik Johansson

******************************************************************************/

#include "partitions.h"

#include "arith.h"
#include "fmprb.h"
#include "math.h"


#define DOUBLE_PREC 53
#define PI 3.141592653589793238462643
#define INV_LOG2 (1.44269504088896340735992468 + 1e-12)
#define HRR_A (1.1143183348516376904 + 1e-12)  /* 44*pi^2/(225*sqrt(3)) */
#define HRR_B (0.0592384391754448833 + 1e-12)  /* pi*sqrt(2)/75 */
#define HRR_C (2.5650996603237281911 + 1e-12)  /* pi*sqrt(2/3) */
#define HRR_D (1.2424533248940001551 + 1e-12)  /* log(2) + log(3)/2 */

static void
_fmpr_sinh(fmpr_t y, const fmpr_t x, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_set_fmpr(t, x);
    fmprb_sinh(t, t, prec);
    fmpr_add(y, fmprb_midref(t), fmprb_radref(t), prec, FMPR_RND_UP);
    fmprb_clear(t);
}


/* Equation (1.8) in the paper */
void
partitions_remainder_bound_fmpr(fmpr_t b, ulong n, ulong N)
{
    fmpr_t A, B, C, t, u;

    fmpr_init(A);
    fmpr_init(B);
    fmpr_init(C);
    fmpr_init(t);
    fmpr_init(u);

    /* bound for 44*pi^2/(225*sqrt(3)) */
    fmpr_set_si_2exp_si(A, 18695160, -24);

    /* bound for pi*sqrt(2)/75 */
    fmpr_set_si_2exp_si(B, 993857, -24);

    /* bound for pi*sqrt(2/3) */
    fmpr_set_si_2exp_si(C, 43035232, -24);

    /* first term: A / sqrt(N) */
    fmpr_sqrt_ui(t, N, FMPRB_RAD_PREC, FMPR_RND_DOWN);
    fmpr_div(b, A, t, FMPRB_RAD_PREC, FMPR_RND_UP);

    /* B * sqrt(N/(n-1)) */
    fmpr_set_ui(t, N);
    fmpr_div_ui(t, t, n - 1, FMPRB_RAD_PREC, FMPR_RND_UP);
    fmpr_sqrt(t, t, FMPRB_RAD_PREC, FMPR_RND_UP);
    fmpr_mul(t, B, t, FMPRB_RAD_PREC, FMPR_RND_UP);

    /* sinh(C*sqrt(n)/N) */
    fmpr_sqrt_ui(u, n, FMPRB_RAD_PREC, FMPR_RND_UP);
    fmpr_div_ui(u, u, N, FMPRB_RAD_PREC, FMPR_RND_UP);
    fmpr_mul(u, C, u, FMPRB_RAD_PREC, FMPR_RND_UP);

    _fmpr_sinh(u, u, FMPRB_RAD_PREC);

    /* second term: B * ... * sinh... */
    fmpr_mul(t, t, u, FMPRB_RAD_PREC, FMPR_RND_UP);
    fmpr_add(b, b, t, FMPRB_RAD_PREC, FMPR_RND_UP);

    fmpr_clear(A);
    fmpr_clear(B);
    fmpr_clear(C);
    fmpr_clear(t);
    fmpr_clear(u);
}


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

static long
partitions_needed_terms(ulong n)
{
    long N;
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
    6469693230UL, 200560490130UL, 7420738134810UL, 304250263527210UL,
    13082761331670030UL, 614889782588491410UL
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


static __inline__ long
log2_ceil(double x)
{
    /* ceil(log2(n)) = bitcount(n-1);
       this is too large if x is a power of two */
    return FLINT_BIT_COUNT((long) x);
}

static long
partitions_prec_bound(ulong n, long k, long N)
{
    long prec;

    prec = partitions_term_bound(n, k);
    prec += log2_ceil(8 * N * (26 * (sqrt(n) / k) + 7 * bound_primes(k) + 22));

    return prec;
}



static void
eval_trig_prod(fmprb_t sum, trig_prod_t prod, long prec)
{
    int i;

    if (prod->prefactor == 0)
    {
        fmprb_zero(sum);
        return;
    }

    if (prec <= DOUBLE_PREC && 0)
    {
    }
    else
    {
        mp_limb_t v;
        fmprb_t t;

        fmprb_init(t);

        fmprb_set_si(sum, prod->prefactor);
        v = n_gcd(FLINT_MAX(prod->sqrt_p, prod->sqrt_q),
                  FLINT_MIN(prod->sqrt_p, prod->sqrt_q));
        prod->sqrt_p /= v;
        prod->sqrt_q /= v;

        if (prod->sqrt_p != 1)
        {
            fmprb_sqrt_ui(t, prod->sqrt_p, prec);
            fmprb_mul(sum, sum, t, prec);
        }

        if (prod->sqrt_q != 1)
        {
            fmprb_sqrt_ui(t, prod->sqrt_q, prec);
            fmprb_div(sum, sum, t, prec);
        }

        for (i = 0; i < prod->n; i++)
        {
            fmpq_t pq;
            *fmpq_numref(pq) = prod->cos_p[i];
            *fmpq_denref(pq) = prod->cos_q[i];
            fmprb_cos_pi_fmpq(t, pq, prec);
            fmprb_mul(sum, sum, t, prec);
        }

        fmprb_clear(t);
    }
}

static void
sinh_cosh_divk_precomp(fmprb_t sh, fmprb_t ch, fmprb_t ex, long k, long prec)
{
    fmprb_t t;
    fmprb_init(t);
    fmprb_set_round(t, ex, prec);
    fmprb_root(ch, t, k, prec);
    /* The second term doesn't need full precision,
       but this doesn't affect performance that much... */
    fmprb_ui_div(t, 1, ch, prec);
    fmprb_sub(sh, ch, t, prec);
    fmprb_add(ch, ch, t, prec);
    fmprb_mul_2exp_si(ch, ch, -1);
    fmprb_mul_2exp_si(sh, sh, -1);
    fmprb_clear(t);
}


void
partitions_hrr_sum_fmprb(fmprb_t x, ulong n, long N0, long N)
{
    trig_prod_t prod;
    fmprb_t acc, C, t1, t2, t3, t4, exp1;
    fmpr_t bound;
    fmpz_t n24;
    long k, prec, res_prec, acc_prec, guard_bits;

    if (n <= 2)
    {
        fmprb_set_ui(x, FLINT_MAX(1, n));
        return;
    }

    /* compute initial precision */
    guard_bits = 2 * FLINT_BIT_COUNT(N) + 32;
    prec = partitions_remainder_bound_log2(n, N0) + guard_bits;
    prec = FLINT_MAX(prec, DOUBLE_PREC);
    res_prec = acc_prec = prec;

    fmprb_init(acc);
    fmprb_init(C);
    fmprb_init(t1);
    fmprb_init(t2);
    fmprb_init(t3);
    fmprb_init(t4);
    fmprb_init(exp1);
    fmpz_init(n24);

    fmprb_zero(x);

    /* n24 = 24n - 1 */
    fmpz_set_ui(n24, n);
    fmpz_mul_ui(n24, n24, 24);
    fmpz_sub_ui(n24, n24, 1);

    /* C = (pi/6) sqrt(24n-1) */
    fmprb_const_pi(t1, prec);
    fmprb_sqrt_fmpz(t2, n24, prec);
    fmprb_mul(t1, t1, t2, prec);
    fmprb_div_ui(C, t1, 6, prec);

    /* exp1 = exp(C) */
    fmprb_exp(exp1, C, prec);

    for (k = N0; k <= N; k++)
    {
        trig_prod_init(prod);
        arith_hrr_expsum_factored(prod, k, n % k);

        if (prod->prefactor != 0)
        {
            if (prec > DOUBLE_PREC)
                prec = partitions_prec_bound(n, k, N);

            /* Compute A_k(n) * sqrt(3/k) * 4 / (24*n-1) */
            prod->prefactor *= 4;
            prod->sqrt_p *= 3;
            prod->sqrt_q *= k;
            eval_trig_prod(t1, prod, prec);
            fmprb_div_fmpz(t1, t1, n24, prec);

            /* Multiply by (cosh(z) - sinh(z)/z) where z = C / k */
            fmprb_set_round(t2, C, prec);
            fmprb_div_ui(t2, t2, k, prec);

            if (k < 35)
                sinh_cosh_divk_precomp(t3, t4, exp1, k, prec);
            else
                fmprb_sinh_cosh(t3, t4, t2, prec);

            fmprb_div(t3, t3, t2, prec);
            fmprb_sub(t2, t4, t3, prec);
            fmprb_mul(t1, t1, t2, prec);

            /* Add to accumulator */
            fmprb_add(acc, acc, t1, acc_prec);

            if (acc_prec > 2 * prec + 32)
            {
                fmprb_add(x, x, acc, res_prec);
                fmprb_zero(acc);
                acc_prec = prec + 32;
            }
        }
    }

    fmprb_add(x, x, acc, res_prec);

    fmpr_init(bound);
    partitions_rademacher_bound(bound, n, N);
    /* printf("addery %lu %lu: ", n, N); fmpr_printd(bound, 20); printf("\n"); */
    fmpr_add(fmprb_radref(x), fmprb_radref(x), bound, FMPRB_RAD_PREC, FMPR_RND_UP);
    fmpr_clear(bound);

    fmpz_clear(n24);
    fmprb_clear(acc);
    fmprb_clear(exp1);
    fmprb_clear(C);
    fmprb_clear(t1);
    fmprb_clear(t2);
    fmprb_clear(t3);
    fmprb_clear(t4);
}

void
partitions_fmpz_ui(fmpz_t p, ulong n)
{
    fmprb_t x;
    fmprb_init(x);
    partitions_hrr_sum_fmprb(x, n, 1, partitions_needed_terms(n));

    if (!fmprb_get_unique_fmpz(p, x))
    {
        printf("not unique!\n");
        fmprb_printd(x, 50);
        printf("\n");
        abort();
    }

    fmprb_clear(x);
}

