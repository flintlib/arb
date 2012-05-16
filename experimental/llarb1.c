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
    long exp;
}
arb1_struct;

typedef arb1_struct arb1_t[1];
typedef arb1_struct * arb1_ptr;
typedef const arb1_struct * arb1_srcptr;


void
arb1_init(arb1_t x, long limbs)
{
    x->d = malloc(sizeof(mp_limb_t) * limbs);
    x->d[0] = 0;
    x->alloc = limbs;
    x->rad = 0UL;
    x->size = 1;
    x->exp = 0;
}

void
arb1_clear(arb1_t x)
{
    free(x->d);
}

void
arb1_randtest(arb1_t x)
{
    mp_limb_t n, m;

    mpn_random(&n, 1);
    n = 1 + (n % x->alloc);

    mpn_random2(x->d, n);
    x->size = n;

    mpn_random(&n, 1);
    if (n % 2)
        x->size *= -1;

    mpn_random2(&n, 1);
    mpn_random2(&m, 1);
    x->rad = (n << 32) | (m >> 32);

    mpn_random(&n, 1);
    x->exp = (n % (2*x->alloc)) - x->alloc;
}

void
arb1_debug(arb1_t x)
{
    printf("alloc=%ld size=%ld exp=%ld rad=%lu", x->alloc, x->size, x->exp, x->rad);
    mpn_debug(x->d, FLINT_ABS(x->size));
    printf("\n");
}

void
arb1_get_rad_mpfr(mpfr_t r, arb1_srcptr x)
{
    mpfr_set_ui_2exp(r, x->rad, FLINT_BITS * (x->exp - FLINT_ABS(x->size)), MPFR_RNDU);
}

long
arb1_get_mpz_2exp(mpz_t a, arb1_srcptr x)
{
    if ((FLINT_ABS(x->size) == 1) && (x->d[0] == 0))
    {
        mpz_set_ui(a, 0);
        return 0;
    }
    else
    {
        mpz_realloc2(a, FLINT_ABS(x->size) * FLINT_BITS);
        mpn_copyi(a->_mp_d, x->d, FLINT_ABS(x->size));
        a->_mp_size = x->size;
        return (x->exp - FLINT_ABS(x->size)) * FLINT_BITS;
    }
}

void
arb1_get_mpfr_interval(mpfr_t a, mpfr_t b, arb1_srcptr x)
{
    long e;
    mpfr_t r;
    mpz_t t;
    mpz_init(t);
    mpfr_init2(r, FLINT_BITS);

    e = arb1_get_mpz_2exp(t, x);

    mpfr_set_z_2exp(a, t, e, MPFR_RNDD);
    mpfr_set_z_2exp(b, t, e, MPFR_RNDU);

    arb1_get_rad_mpfr(r, x);

    mpfr_sub(a, a, r, MPFR_RNDD);
    mpfr_add(b, b, r, MPFR_RNDU);

    mpfr_clear(r);
    mpz_clear(t);
}

int
arb1_contains_mpfr(arb1_srcptr x, const mpfr_t f)
{
    mpfr_t a, b;
    int result;

    mpfr_init2(a, FLINT_BITS * FLINT_ABS(x->size));
    mpfr_init2(b, FLINT_BITS * FLINT_ABS(x->size));

    arb1_get_mpfr_interval(a, b, x);

    result = (mpfr_cmp(a, f) <= 0) && (mpfr_cmp(f, b) <= 0);

    mpfr_clear(a);
    mpfr_clear(b);

    return result;
}


#define TOP_LIMBS(xx,xxn,nnn) ((xx)->d + (xxn) - (nnn))

void
arb1_mul(arb1_ptr z, arb1_srcptr x, arb1_srcptr y)
{
    long zn, xn, yn, zexp, len, topzero, shift, radexp;
    mp_limb_t xrad, yrad, zrad, t1, t0, u1, u0, xtop, ytop;
    mp_limb_t tmp[64];

    zn = z->alloc;

    xn = FLINT_ABS(x->size);
    yn = FLINT_ABS(y->size);

    zn = FLINT_MIN(zn, xn);
    zn = FLINT_MIN(zn, yn);

    xrad = (xn == zn) ? x->rad : 2UL;
    yrad = (yn == zn) ? y->rad : 2UL;

    xtop = TOP_LIMBS(x, xn, 1)[0];
    ytop = TOP_LIMBS(y, yn, 1)[0];

    if (zn == 1)
    {
        mp_limb_t u1, u0, a2, a1, a0, b1, b0;

        if (xtop != ((mp_limb_t) -1))
        {
            umul_ppmm(a1, a0, xtop + 1, yrad);
        }
        else
        {
            a1 = yrad;
            a0 = 0;
        }

        if (ytop != ((mp_limb_t) -1))
        {
            umul_ppmm(b1, b0, ytop + 1, xrad);
        }
        else
        {
            b1 = xrad;
            b0 = 0;
        }

        add_sssaaaaaa(a2, a1, a0, 0, a1, a0, 0, b1, b0);

        umul_ppmm(b1, b0, xrad, yrad);
        add_sssaaaaaa(a2, a1, a0, a2, a1, a0, 0, b1, b0);

        umul_ppmm(u1, u0, xtop, ytop);

        if (a2 != 0 || (a1 >= ((mp_limb_t) -2)))
        {
            z->rad = a2 + 1;
            z->d[0] = 0;
            z->size = 1;
            z->exp = x->exp + y->exp + 2;
        }
        else if (a1 != 0 || (a0 >= ((mp_limb_t) -2)))
        {
            z->rad = a1 + 2;
            z->d[0] = u1;
            z->size = ((x->size ^ y->size) < 0) ? -1 : 1;
            z->exp = x->exp + y->exp;
        }
        else
        {
            if (u1 == 0)
            {
                z->rad = a0;
                z->d[0] = u0;
                z->size = ((x->size ^ y->size) < 0) ? -1 : 1;
                z->exp = x->exp + y->exp - 1;
            }
            else
            {
                z->rad = 2;
                z->d[0] = u1;
                z->size = ((x->size ^ y->size) < 0) ? -1 : 1;
                z->exp = x->exp + y->exp;
            }
        }

        return;
    }

    /* shift down to allow computing xtop*yrad + ytop*xrad + roundoff
       without overflowing two limbs */
#if 1
    xtop = (xtop >> 2) + 2;
    ytop = (ytop >> 2) + 2;
#else
    xrad = (xrad >> 2) + 2;
    yrad = (yrad >> 2) + 2;
#endif

    umul_ppmm(t1, t0, xrad, ytop);
    umul_ppmm(u1, u0, yrad, xtop);
    add_ssaaaa(t1, t0, t1, t0, u1, u0);

#define MAX_UNCORRECTED_RAD ((1UL << (FLINT_BITS - 2)) - 2)

    radexp = 0;

    if (t1 >= MAX_UNCORRECTED_RAD)
    {
        zrad = 4;
        radexp = zn + 1;
    }
    else if (t1 != 0 || t0 >= MAX_UNCORRECTED_RAD)
    {
        zrad = t1 * 4 + 6;
        radexp = zn;
    }
    else
    {
        zrad = t0 * 4 + 2;
        radexp = zn - 1;
    }

    mpn_mul_n(tmp, TOP_LIMBS(x, xn, zn), TOP_LIMBS(y, yn, zn), zn);

    topzero = (tmp[2 * zn - 1] == 0);
    len = 2*zn - topzero - radexp;

    if (len > z->alloc)
    {
        radexp++;
        zrad = 2UL;
        len = 2*zn - topzero - radexp;
    }

    mpn_copyi(z->d, tmp + radexp, 2*zn - topzero - radexp);

    z->size = ((x->size ^ y->size) < 0) ? -len : len;
    z->rad = zrad;
    z->exp = x->exp + y->exp - topzero;
}

void
arb1_add(arb1_ptr z, arb1_srcptr x, arb1_srcptr y)
{
    long zn, xn, yn, xexp, yexp, distance;
    int sign, subtraction;
    mp_limb_t cy, xrad, yrad, rb, ra;

    xexp = x->exp;
    yexp = y->exp;

    if (xexp < yexp)
    {
        arb1_srcptr __t;
        long __u;
        __t = x; x = y; y = __t;
        __u = xexp; xexp = yexp; yexp = __u;
    }

    xn = x->size;
    yn = y->size;
    zn = z->alloc;

    subtraction = (xn ^ yn) < 0;
    sign = xn < 0;

    xn = FLINT_ABS(xn);
    yn = FLINT_ABS(yn);
    zn = FLINT_MIN(zn, xn);
    distance = xexp - yexp;

    /* same signs (addition) */
    if (!subtraction)
    {
        if (distance < zn)
        {   /* y overlaps with the result */
            zn = FLINT_MIN(zn, distance + yn);
            cy = mpn_add(z->d, TOP_LIMBS(x, xn, zn), zn,
                               TOP_LIMBS(y, yn, zn - distance), zn - distance);
            yrad = (zn == (distance + yn)) ? y->rad : 2UL;
        }
        else
        {   /* y does not overlap with result */
            if (x != y)
                mpn_copyi(z->d, TOP_LIMBS(x, xn, zn), zn);
            cy = 0;
            yrad = 2UL;
        }

        xrad = (zn == xn) ? x->rad : 2UL;
        add_ssaaaa(rb, ra, 0, xrad, 0, yrad);

        z->size = zn;
        z->rad = ra;
        z->exp = xexp;

        /* we need to remove one limb */
        if (rb != 0 || (cy != 0 && z->alloc <= zn))
        {
            if (zn != 1)
                mpn_copyi(z->d, z->d + 1, zn - 1);
            z->d[zn - 1] = cy;
            z->size = zn - 1 + cy;
            z->rad = rb + 2UL;
            z->exp += (cy != 0);

            if (z->size == 0)
            {
                z->size = 1;
                z->exp++;
            }
        }
        else if (cy != 0)  /* carry, but it fits */
        {
            z->d[zn] = cy;
            z->size++;
            z->exp++;
        }
    }
    else
    {   /* different signs (subtraction) */

        if (distance == 0)
        {   /* y overlaps with the result; both are aligned */
            mp_ptr xp, yp;

            zn = FLINT_MIN(zn, yn);
            xp = TOP_LIMBS(x, xn, zn);
            yp = TOP_LIMBS(y, yn, zn);

            if (mpn_cmp(xp, yp, zn) >= 0)
            {
                mpn_sub_n(z->d, xp, yp, zn);
            }
            else
            {
                mpn_sub_n(z->d, yp, xp, zn);
                sign = !sign;
            }

            yrad = (zn == yn) ? y->rad : 2UL;
        }
        else if (distance < zn)
        {   /* y overlaps with the result */
            zn = FLINT_MIN(zn, distance + yn);
            mpn_sub(z->d, TOP_LIMBS(x, xn, zn), zn,
                          TOP_LIMBS(y, yn, zn - distance), zn - distance);
            yrad = (zn == (distance + yn)) ? y->rad : 2UL;
        }
        else
        {   /* y does not overlap with result */
            if (x != y)
                mpn_copyi(z->d, TOP_LIMBS(x, xn, zn), zn);
            yrad = 2UL;
        }

        xrad = (zn == xn) ? x->rad : 2UL;
        add_ssaaaa(rb, ra, 0, xrad, 0, yrad);

        z->rad = ra;

        while (zn > 1 && z->d[zn - 1] == 0)
        {
            zn--;
            xexp--;
        }

        z->exp = xexp;
        z->size = zn;

        if (rb != 0)
        {
            if (zn != 1)
                mpn_copyi(z->d, z->d + 1, zn - 1);
            z->rad = rb + 2UL;
            z->size = zn - 1;
        }

        if (z->size == 0)
        {
            z->d[0] = 0;
            z->size = 1;
            z->exp++;
        }
    }

    if (sign)
        z->size = -z->size;
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
        arb1_t x, y, z;
        mpfr_t t, u, v;

        arb1_init(x, 1 + n_randint(state, 10));
        arb1_init(y, 1 + n_randint(state, 10));
        arb1_init(z, 1 + n_randint(state, 10));

        arb1_randtest(x);
        arb1_randtest(y);
        arb1_randtest(z);

        mpfr_init2(t, 2000);
        mpfr_init2(u, 2000);
        mpfr_init2(v, 2000);

        if (n_randint(state, 2))
            arb1_get_mpfr_interval(t, v, x);
        else
            arb1_get_mpfr_interval(v, t, x);

        if (n_randint(state, 2))
            arb1_get_mpfr_interval(u, v, y);
        else
            arb1_get_mpfr_interval(v, u, y);

        if (FLINT_ABS(x->size) >= 1 && FLINT_ABS(y->size) >= 1 && z->alloc >= 1)
        {
            arb1_mul(z, x, y);
            mpfr_mul(v, t, u, MPFR_RNDN);

            if (!arb1_contains_mpfr(z, v))
            {
                printf("FAIL! [%ld]\n", iter);
                arb1_debug(x);
                arb1_debug(y);
                arb1_debug(z);

                mpfr_debug(t);
                mpfr_debug(u);
                mpfr_debug(v);

                abort();
            }
        }

        mpfr_clear(t);
        mpfr_clear(u);
        mpfr_clear(v);

        arb1_clear(x);
        arb1_clear(y);
        arb1_clear(z);
    }

    flint_randclear(state);
    return EXIT_SUCCESS;
}

int testadd()
{
    long iter;
    flint_rand_t state;

    printf("testadd....");
    fflush(stdout);

    flint_randinit(state);
    _flint_rand_init_gmp(state);

    for (iter = 0; iter < 1000000; iter++)
    {
        arb1_t x, y, z;
        mpfr_t t, u, v;

        arb1_init(x, 1 + n_randint(state, 10));
        arb1_init(y, 1 + n_randint(state, 10));
        arb1_init(z, 1 + n_randint(state, 10));

        arb1_randtest(x);
        arb1_randtest(y);
        arb1_randtest(z);

        mpfr_init2(t, 2000);
        mpfr_init2(u, 2000);
        mpfr_init2(v, 2000);

        if (n_randint(state, 2))
            arb1_get_mpfr_interval(t, v, x);
        else
            arb1_get_mpfr_interval(v, t, x);

        if (n_randint(state, 2))
            arb1_get_mpfr_interval(u, v, y);
        else
            arb1_get_mpfr_interval(v, u, y);

        if (1)
        {
            arb1_add(z, x, y);
            mpfr_add(v, t, u, MPFR_RNDN);

            if (!arb1_contains_mpfr(z, v))
            {
                printf("FAIL! [%ld]\n", iter);
                arb1_debug(x);
                arb1_debug(y);
                arb1_debug(z);

                mpfr_debug(t);
                mpfr_debug(u);
                mpfr_debug(v);

                abort();
            }
        }

        mpfr_clear(t);
        mpfr_clear(u);
        mpfr_clear(v);

        arb1_clear(x);
        arb1_clear(y);
        arb1_clear(z);
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
    arb1_t x[N], y[N], z[N];
    long i, j, k, limbs;

    flint_rand_t state;
    flint_randinit(state);

    printf("\n\n\n");

    for (limbs = 1; limbs <= 15; limbs++)
    {
        for (i = 0; i < N; i++)
        {
            arb1_init(x[i], limbs);
            arb1_init(y[i], limbs);
            arb1_init(z[i], limbs);

            do { arb1_randtest(x[i]); } while (x[i]->size != limbs);
            do { arb1_randtest(y[i]); } while (y[i]->size != limbs);
            do { arb1_randtest(z[i]); } while (z[i]->size != limbs);
        }

        printf("limbs = %ld   ", limbs);

        TIMEIT_START
        for (i = 0; i < N; i++)
            arb1_add(z[i], x[i], y[i]);
        TIMEIT_STOP
        printf("    ");

        TIMEIT_START
        for (i = 0; i < N; i++)
            arb1_mul(z[i], x[i], y[i]);
        TIMEIT_STOP
        printf("\n");

        for (i = 0; i < N; i++)
        {
            arb1_clear(x[i]);
            arb1_clear(y[i]);
            arb1_clear(z[i]);
        }
    }

    printf("\n\n");

    testmul();
    testadd();

    _fmpz_cleanup();
}
