#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <mpir.h>
#include <mpfr.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "mpn_extras.h"
#include "mpn_inlines.h"
#include "ufloat.h"

#define _MPR_BITS_TO_LIMBS(b) (((b) + FLINT_BITS - 1) / FLINT_BITS)

/* TODO: add target precision field? */

typedef struct
{
    mp_ptr d;
    long exp;
    long size;
    long alloc;
    int sign;
}
mpr_struct;

typedef mpr_struct mpr_t[1];
typedef mpr_struct * mpr_ptr;
typedef const mpr_struct * mpr_srcptr;

#define MPR_SIGN_PLUS 1
#define MPR_SIGN_MINUS -1



/* Memory management *********************************************************/

static __inline__ void
mpr_init(mpr_t x)
{
    x->d = NULL;
    x->exp = 0;
    x->sign = MPR_SIGN_PLUS;
    x->alloc = 0;
    x->size = 0;
}

static __inline__ void
mpr_clear(mpr_t x)
{
    free(x->d);
}

static __inline__ void
mpr_fit_limbs(mpr_t x, mp_size_t n)
{
    if (n > x->alloc)
    {
        if (x->alloc == 0)
            x->d = (mp_ptr) malloc(n * sizeof(mp_limb_t));
        else
            x->d = (mp_ptr) realloc(x->d, n * sizeof(mp_limb_t));
        x->alloc = n;
    }
}

static __inline__ void
mpr_fit_bits(mpr_t x, mp_bitcnt_t bits)
{
    mpr_fit_limbs(x, _MPR_BITS_TO_LIMBS(bits));
}

void mpr_normalize(mpr_t x);

int mpr_is_normalized(mpr_t x);

static __inline__ void
mpr_zero(mpr_t x)
{
    x->size = 0;
    x->exp = 0;
    x->sign = MPR_SIGN_PLUS;
}

static __inline__ int
mpr_is_zero(const mpr_t x)
{
    return x->size == 0 && x->exp == 0;
}

static __inline__ void
mpr_set(mpr_t y, const mpr_t x)
{
    if (y != x)
    {
        mpr_fit_limbs(y, x->size);
        mpn_copyi(y->d, x->d, x->size);

        y->exp = x->exp;
        y->size = x->size;
        y->sign = x->sign;
    }
}

static __inline__ void
mpr_swap(mpr_t y, mpr_t x)
{
    if (y != x)
    {
        mpr_struct t = *x;
        *x = *y;
        *y = t;
    }
}

static __inline__ void
mpr_abs(mpr_t y, const mpr_t x)
{
    if (y != x)
        mpr_set(y, x);

    if (y->sign == MPR_SIGN_MINUS)
        y->sign = MPR_SIGN_PLUS;
}

static __inline__ void
mpr_neg(mpr_t y, const mpr_t x)
{
    if (y != x)
        mpr_set(y, x);

    if (!mpr_is_zero(y))
        y->sign = -y->sign;
}


/* Conversions ***************************************************/

static __inline__ void mpr_get_mpfr_wrapper(mpfr_t s, mpr_srcptr x)
{
    s->_mpfr_d = (mp_ptr) x->d;
    s->_mpfr_exp = x->exp;
    s->_mpfr_prec = x->size * FLINT_BITS;
    s->_mpfr_sign = (x->sign == MPR_SIGN_PLUS) ? 1 : -1;

    if (mpr_is_zero(x))
    {
        s->_mpfr_exp = __MPFR_EXP_ZERO;
    }
}


long _mpr_set_mpfr(mp_ptr y, const mpfr_t x, mp_size_t n);
void _mpr_get_mpfr(mpfr_t y, mp_srcptr x, long exp, mp_size_t n, mpfr_rnd_t rnd);
void _mpr_get_mpfr_signed(mpfr_t y, mp_srcptr x, long exp,
    mp_size_t n, int sign, mpfr_rnd_t rnd);

void mpr_set_mpfr(mpr_t y, const mpfr_t x);
void mpr_get_mpfr(mpfr_t y, const mpr_t x, mpfr_rnd_t rnd);

void mpr_set_ui(mpr_t y, ulong x);
void mpr_set_si(mpr_t y, long x);
void mpr_set_mpz(mpr_t y, const mpz_t x);
void mpr_set_fmpz(mpr_t y, const fmpz_t x);

static __inline__ void
mpr_mul_2exp(mpr_t y, const mpr_t x, long exp)
{
    if (x != y)
        abort();

    if (y->size != 0)
        y->exp += exp;
}

static __inline__ void
mpr_set_ui_2exp(mpr_t y, ulong x, long exp)
{
    mpr_set_ui(y, x);
    mpr_mul_2exp(y, y, exp);
}

static __inline__ void
mpr_set_si_2exp(mpr_t y, long x, long exp)
{
    mpr_set_si(y, x);
    mpr_mul_2exp(y, y, exp);
}

static __inline__ void
mpr_set_mpz_2exp(mpr_t y, const mpz_t x, long exp)
{
    mpr_set_mpz(y, x);
    mpr_mul_2exp(y, y, exp);
}

static __inline__ void
mpr_set_fmpz_2exp(mpr_t y, const fmpz_t x, long exp)
{
    mpr_set_fmpz(y, x);
    mpr_mul_2exp(y, y, exp);
}

void _mpr_randtest(mp_ptr x, flint_rand_t state, mp_size_t n);
void mpr_randtest(mpr_t x, flint_rand_t state, mp_size_t n);

void _mpr_printr(mp_ptr x, mp_size_t n);

static __inline__ void
mpr_debug(const mpr_t x)
{
    long i;

    printf("mpr(sign=%d, exp=%lu, size=%lu, d=[", x->sign, x->exp, x->size);

    for (i = 0; i < x->size; i++)
    {
        printf("%lu", x->d[i]);
        if (i < x->size - 1)
            printf(", ");
    }

    printf("])\n");
}


/* XXX: assumes nonzero! */
static __inline__ void
_mpr_get_ufloat(ufloat_t u, mp_srcptr d, long size, long exp)
{
    /* since we truncate, we need to round up */
    u->man = (d[size - 1] >> (FLINT_BITS - UFLOAT_PREC)) + 1UL;
    u->exp = exp;

    /* adjust for carry (very unlikely!) */
    if (u->man >= (1UL << UFLOAT_PREC))
    {
        u->man = (u->man >> 1) + 1UL;
        u->exp++;
    }
}

static __inline__ void
mpr_get_ufloat(ufloat_t u, const mpr_t x)
{
    _mpr_get_ufloat(u, x->d, x->size, x->exp);
}

static __inline__ int
mpr_equal(const mpr_t x, const mpr_t y)
{
    if (x->sign == y->sign && x->exp == y->exp && x->size == y->size)
    {
        return mpn_cmp(x->d, y->d, x->size) == 0;
    }

    return 0;
}

/* Exact arithmetic **********************************************************/

void mpr_mul_exact(mpr_ptr z, mpr_srcptr x, mpr_srcptr y);

/* Arithmetic using MPFR *****************************************************/

long mpr_add_using_mpfr(mpr_ptr z, mpr_srcptr x, mpr_srcptr y, mp_bitcnt_t prec, mpfr_rnd_t rnd);
long mpr_sub_using_mpfr(mpr_ptr z, mpr_srcptr x, mpr_srcptr y, mp_bitcnt_t prec, mpfr_rnd_t rnd);
long mpr_mul_using_mpfr(mpr_ptr z, mpr_srcptr x, mpr_srcptr y, mp_bitcnt_t prec, mpfr_rnd_t rnd);

long mpr_mul_ui_using_mpfr(mpr_ptr z, mpr_srcptr x, ulong y, mp_bitcnt_t prec, mpfr_rnd_t rnd);
long mpr_mul_si_using_mpfr(mpr_ptr z, mpr_srcptr x, long y, mp_bitcnt_t prec, mpfr_rnd_t rnd);

long mpr_addmul_using_mpfr(mpr_ptr z, mpr_srcptr x, mpr_srcptr y, mp_bitcnt_t prec, mpfr_rnd_t rnd);
long mpr_submul_using_mpfr(mpr_ptr z, mpr_srcptr x, mpr_srcptr y, mp_bitcnt_t prec, mpfr_rnd_t rnd);
long mpr_addmul_ui_using_mpfr(mpr_ptr z, mpr_srcptr x, ulong y, mp_bitcnt_t prec, mpfr_rnd_t rnd);
long mpr_submul_ui_using_mpfr(mpr_ptr z, mpr_srcptr x, ulong y, mp_bitcnt_t prec, mpfr_rnd_t rnd);



/* Experimental arithmetic ***************************************************/

int _mpr_addd_n(mp_ptr z, mp_srcptr x, mp_srcptr y, mp_bitcnt_t shift, mp_size_t n);
int _mpr_addd_using_mpfr(mp_ptr z, mp_size_t zn, mp_srcptr x, mp_size_t xn, mp_srcptr y, mp_size_t yn, long shift);
int _mpr_muld_n(mp_ptr z, mp_srcptr x, mp_srcptr y, mp_size_t n);
int _mpr_muld_using_mpfr(mp_ptr z, mp_size_t zn, mp_srcptr x, mp_size_t xn, mp_srcptr y, mp_size_t yn);
