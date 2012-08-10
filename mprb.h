#include "mpr.h"
#include "ufloat.h"

typedef struct
{
    mpr_struct mid;
    ufloat_struct rad;
}
mprb_struct;

typedef mprb_struct mprb_t[1];
typedef mprb_struct * mprb_ptr;
typedef const mprb_struct * mprb_srcptr;

#define MPRB_SIGN_PLUS 0
#define MPRB_SIGN_MINUS 1

void mprb_init(mprb_t x, long bits);
void mprb_clear(mprb_t x);

void mprb_debug(const mprb_t x);
void mprb_randtest(mprb_t x, flint_rand_t state, long emin, long emax);

long mprb_get_mid_mpz_2exp(mpz_t a, const mprb_t x);

void mprb_set_mpfr(mprb_t x, const mpfr_t v);
void mprb_get_rad_mpfr(mpfr_t r, const mprb_t x);
void mprb_get_mid_mpfr(mpfr_t x, const mprb_t v, mpfr_rnd_t rnd);
void mprb_get_interval_mpfr(mpfr_t a, mpfr_t b, const mprb_t x);
int mprb_contains_mpfr(const mprb_t x, const mpfr_t f);




void mprb_get_lower_bound_ufloat(ufloat_t u, const mprb_t x);

static __inline__ int
mprb_is_exact(const mprb_t x)
{
    return ufloat_is_zero(&x->rad);
}

int mprb_contains_zero(const mprb_t x);


static __inline__ long
mprb_ulp_exp(const mprb_t x)
{
    return x->mid.exp - x->mid.size * FLINT_BITS;
}

void mprb_add(mprb_t z, const mprb_t x, const mprb_t y);
void mprb_mul(mprb_t z, const mprb_t x, const mprb_t y);


/* Add r*2^exp to the radius of x */
static __inline__ void
mprb_add_rad_ui_2exp(mprb_t x, mp_limb_t r, long exp)
{
    ufloat_t t;
    ufloat_set_ui_2exp(t, r, exp);
    ufloat_add(&x->rad, &x->rad, t);
}

/* Add 2^exp to the radius of x */
static __inline__ void
mprb_add_rad_2exp(mprb_t x, long exp)
{
    ufloat_add_2exp(&x->rad, &x->rad, exp);
}
