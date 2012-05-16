#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "mpn_extras.h"
#include "profiler.h"

typedef struct
{
    mp_ptr d;
    mp_limb_t rad;
    long size;
    long alloc;
    long sign;
    long exp;
}
arb2_struct;

typedef arb2_struct arb2_t[1];
typedef arb2_struct * arb2_ptr;
typedef const arb2_struct * arb2_srcptr;

#define arb2_SIGN_PLUS   0
#define arb2_SIGN_MINUS  1


static __inline__ void
arb2_init(arb2_t x, long bits)
{
    long limbs = (bits + FLINT_BITS - 1) / FLINT_BITS;

    x->d = calloc(limbs, sizeof(mp_limb_t));
    x->exp = 0;
    x->sign = arb2_SIGN_PLUS;
    x->alloc = limbs;
    x->size = limbs;
    x->rad = 0UL;
}

static __inline__ void
arb2_clear(arb2_t x)
{
    free(x->d);
}

void
arb2_set_mpfr(arb2_t x, const mpfr_t v)
{
    long xprec, vprec, prec;

    xprec = x->alloc;
    vprec = (v->_mpfr_prec + FLINT_BITS - 1) / FLINT_BITS;
    prec = FLINT_MIN(xprec, vprec);

    if (xprec >= vprec)
    {
        long s = xprec - vprec;
        mpn_copyi(x->d + s, v->_mpfr_d, prec);
        if (s != 0)
            mpn_zero(x->d, s);
        x->rad = 0UL;
    }
    else
    {
        mpn_copyi(x->d, v->_mpfr_d + vprec - xprec, prec);
        x->rad = 1UL;
    }

    x->exp = v->_mpfr_exp;
    x->sign = (v->_mpfr_sign == 1) ? 0 : 1;
    x->size = x->alloc;
}

void
arb2_get_rad_mpfr(mpfr_t r, const arb2_t x)
{
    mpfr_set_ui_2exp(r, x->rad, x->exp - FLINT_BITS * x->size, MPFR_RNDU);
}

/* todo: handle 0 */
long
arb2_get_mpz_2exp(mpz_t a, const arb2_t x)
{
    if ((x->size == 1) && (x->d[0] == 0))
        mpz_set_ui(a, 0);
    else
    {
        mpz_realloc2(a, x->size * FLINT_BITS);
        mpn_copyi(a->_mp_d, x->d, x->size);
        a->_mp_size = (x->sign == arb2_SIGN_PLUS) ? x->size : -(x->size);
        return x->exp - x->size * FLINT_BITS;
    }
}

void
arb2_get_mpfr_interval(mpfr_t a, mpfr_t b, const arb2_t x)
{
    long e;
    mpfr_t r;
    mpz_t t;
    mpz_init(t);
    mpfr_init2(r, FLINT_BITS);

    e = arb2_get_mpz_2exp(t, x);
    mpfr_set_z_2exp(a, t, e, MPFR_RNDD);
    mpfr_set_z_2exp(b, t, e, MPFR_RNDU);

    arb2_get_rad_mpfr(r, x);
    mpfr_sub(a, a, r, MPFR_RNDD);
    mpfr_add(b, b, r, MPFR_RNDU);

    mpfr_clear(r);
    mpz_clear(t);
}

int
arb2_contains_mpfr(const arb2_t x, const mpfr_t f)
{
    mpfr_t a, b;
    int result;

    mpfr_init2(a, FLINT_BITS * x->size);
    mpfr_init2(b, FLINT_BITS * x->size);

    arb2_get_mpfr_interval(a, b, x);

    result = (mpfr_cmp(a, f) <= 0) && (mpfr_cmp(f, b) <= 0);

    mpfr_clear(a);
    mpfr_clear(b);

    return result;
}

void
arb2_randtest(arb2_t x, flint_rand_t state)
{
    long n;

    n = n_randint(state, x->alloc) + 1;

    if (n_randint(state, 2))
        mpn_random2(x->d, n);
    else
        mpn_random(x->d, n);

    x->d[n - 1] |= (1UL << (FLINT_BITS - 1));

    /* unnormalised?
    if (n == 1 && n_randint(state, 50) == 0)
        x->d[n - 1] = n_randtest(state);
    */

    x->rad = n_randtest(state);
    x->size = n;

    x->exp = (int) n_randint(state, 100) - 50;
}

void
arb2_debug(const arb2_t x)
{
    long e;
    mpz_t t;

    mpz_init(t);
    e = arb2_get_mpz_2exp(t, x);

    gmp_printf("{mid=%Zx, exp=%ld, rad=%lu, size=%ld, alloc=%ld}\n",
        t, e, x->rad, (long) x->size, (long) x->alloc);

    mpz_clear(t);
}


#define RAD_MUL_SHIFT 8
#define SLOPPY_ERROR 0


void arb2_mul(arb2_t z, const arb2_t x, const arb2_t y)
{
    mp_limb_t xrad, yrad;
    mp_srcptr xptr, yptr;
    long size;

    size = FLINT_MIN(x->size, y->size);
    size = FLINT_MIN(size, z->alloc);
    xrad = (size == x->size) ? x->rad : 2UL;
    yrad = (size == y->size) ? y->rad : 2UL;
    xptr = x->d + x->size - size;
    yptr = y->d + y->size - size;

    if (size == 1)
    {
        mp_limb_t u1, u0, a2, a1, a0, b1, b0;

        umul_ppmm(a1, a0, xptr[0], yrad);
        umul_ppmm(b1, b0, yptr[0], xrad);
        add_sssaaaaaa(a2, a1, a0, 0, a1, a0, 0, b1, b0);

        umul_ppmm(b1, b0, xrad, yrad);
        add_sssaaaaaa(a2, a1, a0, a2, a1, a0, 0, b1, b0);

        umul_ppmm(u1, u0, xptr[0], yptr[0]);

        if ((u1 >> (FLINT_BITS - 1)) == 0)
        {
            u1 = (u1 << 1) | (u0 >> (FLINT_BITS - 1));
            u0 = (u0 << 1);
            a2 = (a2 << 1) | (a1 >> (FLINT_BITS - 1));
            a1 = (a1 << 1) | (a0 >> (FLINT_BITS - 1));
            a0 = (a0 << 1);
            z->exp = x->exp + y->exp - 1;
        }
        else
        {
            z->exp = x->exp + y->exp;
        }

        /* round up */
        add_ssaaaa(a2, a1, a2, a1, 0, (u0 != 0) + (a0 != 0));

        if (a2 == 0)
        {
            z->rad = a1;
            z->d[0] = u1;
        }
        else
        {
            /* TODO: should we always return a nonzero midpoint?
               (this requires embiggening the radius appropriately) */
            z->rad = a2 + 2;
            z->d[0] = 0;
            z->exp += FLINT_BITS;
        }

        z->size = 1;
        z->sign = x->sign ^ y->sign;
    }
    else
    {
        mp_limb_t tmp[32];
        mp_limb_t t0, t1, u0, u1, rad;
        mp_limb_t xtop, ytop;
        int top_bit_set, need_limb_shift;

#if SLOPPY_ERROR
        /* this is very fast, but generically doubles the error */
        /* x*yrad + y*xrad < xrad + yrad */
        add_ssaaaa(t1, t0, 0, xrad, 0, yrad);

        /* xrad * yrad < 1 ulp
           + 1 ulp error in the main multiplication */
        add_ssaaaa(t1, t0, t1, t0, 0, 2UL);
#else
        /* upper bound for x*yrad + y*xrad */
        xtop = (xptr[size-1] >> RAD_MUL_SHIFT) + 1UL;
        ytop = (yptr[size-1] >> RAD_MUL_SHIFT) + 1UL;
        umul_ppmm(t1, t0, xtop, y->rad);
        umul_ppmm(u1, u0, ytop, x->rad);
        add_ssaaaa(t1, t0, t1, t0, u1, u0);

        /* scale back radius to the true magnitude */
        t0 = (t0 >> (FLINT_BITS - RAD_MUL_SHIFT)) | (t1 << RAD_MUL_SHIFT);
        t1 = (t1 >> (FLINT_BITS - RAD_MUL_SHIFT));

        /* xrad * yrad < 1 ulp
           + 1 ulp error in the main multiplication
           + 1 ulp error from the shift above */
        add_ssaaaa(t1, t0, t1, t0, 0, 3UL);
#endif

        /* multiply midpoints */
        mpn_mul_n(tmp, xptr, yptr, size);

        top_bit_set = tmp[2 * size - 1] >> (FLINT_BITS-1);

        if (!top_bit_set)
        {
            t1 = (t1 << 1) | (t0 >> (FLINT_BITS - 1));
            t0 = t0 << 1;
        }

        /* remove one limb */
        if (t1)
        {
            need_limb_shift = 1;

            /* +1 ulp error for the rounding of t1 */
            /* +1 ulp error in the main product, which doubles
               if we perform a shift */
            rad = t1 + 2 + !top_bit_set;
        }
        else
        {
            need_limb_shift = 0;
            rad = t0;
        }

        /* copy result to destination */
        if (top_bit_set)
        {
            mpn_copyi(z->d, tmp + size + need_limb_shift, size - need_limb_shift);
        }
        else
        {
            mpn_lshift(z->d, tmp + size + need_limb_shift,
                size - need_limb_shift, 1);
            z->d[0] |= tmp[size - need_limb_shift] >> (FLINT_BITS - 1);
        }

        z->rad = rad;
        z->exp = x->exp + y->exp + (top_bit_set - 1);
        z->size = size - need_limb_shift;
        z->sign = x->sign ^ y->sign;
    }
}


int testit()
{
    long iter;
    flint_rand_t state;

    printf("testit....");
    fflush(stdout);

    flint_randinit(state);

    _flint_rand_init_gmp(state);

    for (iter = 0; iter < 10000; iter++)
    {
        arb2_t x;
        mpfr_t t;

        arb2_init(x, 1 + n_randint(state, 400));
        arb2_randtest(x, state);

        mpfr_init2(t, 2 + n_randint(state, 400));

        mpfr_urandom(t, state->gmp_state, MPFR_RNDN);
        if (n_randint(state, 2))
            mpfr_neg(t, t, MPFR_RNDN);

        mpfr_mul_2exp(t, t, n_randint(state, 20), MPFR_RNDN);
        mpfr_div_2exp(t, t, n_randint(state, 20), MPFR_RNDN);

        arb2_set_mpfr(x, t);

        if (!arb2_contains_mpfr(x, t))
        {
            printf("FAIL!\n");
            arb2_debug(x);
            mpfr_printf("\n%Rf\n", t);
            mpn_debug(t->_mpfr_d, (t->_mpfr_prec + FLINT_BITS - 1) / FLINT_BITS);
            abort();
        }

        mpfr_clear(t);
        arb2_clear(x);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    return EXIT_SUCCESS;
}

void
mpfr_randtest(mpfr_t t, flint_rand_t state)
{
    mpfr_urandom(t, state->gmp_state, MPFR_RNDN);
    if (n_randint(state, 2))
        mpfr_neg(t, t, MPFR_RNDN);
    mpfr_mul_2exp(t, t, n_randint(state, 50), MPFR_RNDN);
    mpfr_div_2exp(t, t, n_randint(state, 50), MPFR_RNDN);
}

void
mpfr_debug(mpfr_t t)
{
    mpfr_printf("\n%.40Rg\n", t);
/*    mpn_debug(t->_mpfr_d, (t->_mpfr_prec + FLINT_BITS - 1) / FLINT_BITS); */
    printf("\n\n");
}

int testmul()
{
    long iter;
    flint_rand_t state;

    printf("testmul....");
    fflush(stdout);

    flint_randinit(state);

    _flint_rand_init_gmp(state);

    for (iter = 0; iter < 1000000; iter++)
    {
        arb2_t x, y, z;
        mpfr_t t, u, v;

        arb2_init(x, 65 + n_randint(state, 300));
        arb2_init(y, 65 + n_randint(state, 300));
        arb2_init(z, 65 + n_randint(state, 300));

        do { arb2_randtest(x, state); } while (x->size == 0);
        do { arb2_randtest(y, state); } while (y->size == 0);
        do { arb2_randtest(z, state); } while (z->size == 0);

/*
        x->rad = -1;
        y->rad = -1;
*/

        mpfr_init2(t, 1000);
        mpfr_init2(u, 1000);
        mpfr_init2(v, 1000);

        if (n_randint(state, 2))
            arb2_get_mpfr_interval(t, v, x);
        else
            arb2_get_mpfr_interval(v, t, x);

        if (n_randint(state, 2))
            arb2_get_mpfr_interval(u, v, y);
        else
            arb2_get_mpfr_interval(v, u, y);

        arb2_mul(z, x, y);
        mpfr_mul(v, t, u, MPFR_RNDN);

        if (!arb2_contains_mpfr(z, v))
        {
            printf("FAIL!\n");
            arb2_debug(x);
            arb2_debug(y);
            arb2_debug(z);

            mpfr_debug(t);
            mpfr_debug(u);
            mpfr_debug(v);

            abort();
        }

        mpfr_clear(t);
        mpfr_clear(u);
        mpfr_clear(v);

        arb2_clear(x);
        arb2_clear(y);
        arb2_clear(z);
    }

    flint_randclear(state);
    return EXIT_SUCCESS;
}


void print_time(timeit_t t, long reps)
{
    double s;
    double ms;
    double us;
    double ns;

    ms = (double) t->cpu / reps;
    us = ms * 1000;
    ns = us * 1000;
    s = ms / 1000;

    if (s > 10)
        printf("%ld s", (long) s);
    else if (ms > 10)
        printf("%ld ms", (long) ms);
    else if (us > 10)
        printf("%ld us", (long) us);
    else
        printf("%ld ns", (long) ns);
}

#define TIMEIT_START                                        \
    {                                                       \
        long reps = 1;                                      \
        while (1)                                           \
        {                                                   \
            timeit_t t0;                                    \
            timeit_start(t0);                               \
            long _k;                                        \
            for (_k = 0; _k < reps; _k++) {                 \

#define TIMEIT_STOP                                         \
            }                                               \
            timeit_stop(t0);                                \
            if (t0->cpu >= 100)                             \
            {                                               \
                print_time(t0, reps); fflush(stdout);       \
                break;                                      \
            }                                               \
            reps *= 10;                                     \
        }                                                   \
    }


int main()
{
#define N 1000
    arb2_t x[N], y[N], z[N];
    long i, j, k, limbs;

    flint_rand_t state;
    flint_randinit(state);

    testmul();

    printf("\n\n\n");

    for (limbs = 1; limbs <= 15; limbs++)
    {
        for (i = 0; i < N; i++)
        {
            arb2_init(x[i], limbs * FLINT_BITS);
            arb2_init(y[i], limbs * FLINT_BITS);
            arb2_init(z[i], limbs * FLINT_BITS);

            do { arb2_randtest(x[i], state); } while (x[i]->size != limbs);
            do { arb2_randtest(y[i], state); } while (y[i]->size != limbs);
            do { arb2_randtest(z[i], state); } while (z[i]->size != limbs);
        }

        printf("limbs = %ld   ", limbs);

        TIMEIT_START
        for (i = 0; i < N; i++)
            arb2_mul(z[i], x[i], y[i]);
        TIMEIT_STOP
        printf("\n");

        for (i = 0; i < N; i++)
        {
            arb2_clear(x[i]);
            arb2_clear(y[i]);
            arb2_clear(z[i]);
        }
    }

    printf("\n\n");
    _fmpz_cleanup();
}
