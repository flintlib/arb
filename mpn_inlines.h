#include <mpir.h>
#include "flint.h"
#include "longlong.h"

/* Multiplication ************************************************************/

mp_limb_t mpn_mul_basecase(mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t);

#define MPN_MUL_11 umul_ppmm

#define MPN_MUL_21(r2, r1, r0, a1, a0, b0) \
    do { \
        mp_limb_t t1; \
        umul_ppmm(r1, r0, a0, b0); \
        umul_ppmm(r2, t1, a1, b0); \
        add_ssaaaa(r2, r1, r2, r1, 0, t1); \
    } while (0)

#define MPN_MUL_22(r3, r2, r1, r0, a1, a0, b1, b0) \
    do { \
        mp_limb_t t1, t2, t3; \
        umul_ppmm(r1, r0, a0, b0); \
        umul_ppmm(r2, t1, a1, b0); \
        add_ssaaaa(r2, r1, r2, r1, 0, t1); \
        umul_ppmm(t1, t2, a0, b1); \
        umul_ppmm(r3, t3, a1, b1); \
        add_ssaaaa(r3, t1, r3, t1, 0, t3); \
        add_ssaaaa(r2, r1, r2, r1, t1, t2); \
        r3 += r2 < t1; \
    } while (0)

#define MPN_MUL_N(d, a, b, n) \
    do { \
        if (n == 1) \
            MPN_MUL_11(d[1], d[0], a[0], b[0]); \
        else if (n == 2) \
            MPN_MUL_22(d[3], d[2], d[1], d[0], a[1], a[0], b[1], b[0]); \
        else \
            mpn_mul_n(d, a, b, n); \
    } while (0)
