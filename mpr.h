#include <stdio.h>
#include <stdlib.h>
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


/* Conversions / debugging ***************************************************/
long _mpr_set_mpfr(mp_ptr y, const mpfr_t x, mp_size_t n);
void _mpr_get_mpfr(mpfr_t y, mp_srcptr x, long exp, mp_size_t n, mpfr_rnd_t rnd);
void _mpr_get_mpfr_signed(mpfr_t y, mp_srcptr x, long exp,
    mp_size_t n, int sign, mpfr_rnd_t rnd);

void _mpr_randtest(mp_ptr x, flint_rand_t state, mp_size_t n);

void _mpr_printr(mp_ptr x, mp_size_t n);


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

/* Multiplication ************************************************************/

#define _MPR_MUL_SHIFT(shift, d, dn, s, sn) \
    do { \
        int i; \
        shift = !((s)[(sn)-1] >> (FLINT_BITS-1)); \
        for (i = (dn) - 1; i >= 0; i--) \
            (d)[i] = ((s)[i+(sn)-(dn)] << shift) \
                   | (((s)[i+(sn)-(dn)-1] >> (FLINT_BITS-1)) & shift); \
    } while (0)

#define _MPR_MULD_N_CONST(shift, z, x, y, n) \
    do { \
        mp_limb_t __tmp[2*n]; \
        MPN_MUL_N(__tmp, x, y, n); \
        _MPR_MUL_SHIFT(shift, z, n, __tmp, 2*n); \
    } while (0)

#define _MPR_MULD_N(shift, z, x, y, n) \
    do { \
        if (n == 1) \
            _MPR_MULD_N_CONST(shift, z, x, y, 1); \
        else if (n == 2) \
            _MPR_MULD_N_CONST(shift, z, x, y, 2); \
        else if (n <= 15) \
        { \
            mp_limb_t __tmp[30]; \
            mpn_mul_basecase(__tmp, x, n, y, n); \
            _MPR_MUL_SHIFT(shift, z, n, __tmp, 2*n); \
        } \
        else \
        { \
            mp_ptr __tmp = malloc(2 * n * sizeof(mp_limb_t)); \
            mpn_mul_n(__tmp, x, y, n); \
            _MPR_MUL_SHIFT(shift, z, n, __tmp, 2*n); \
            free(__tmp); \
        } \
    } while (0)

int
_mpr_muld_n(mp_ptr z, mp_srcptr x, mp_srcptr y, mp_size_t n);

int _mpr_muld_using_mpfr(mp_ptr z, mp_size_t zn, mp_srcptr x, mp_size_t xn, mp_srcptr y, mp_size_t yn);

/* Addition ******************************************************************/

int _mpr_addd_n(mp_ptr z, mp_srcptr x, mp_srcptr y, mp_bitcnt_t shift, mp_size_t n);

int _mpr_addd_using_mpfr(mp_ptr z, mp_size_t zn, mp_srcptr x, mp_size_t xn, mp_srcptr y, mp_size_t yn, long shift);
