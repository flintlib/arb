#include "mpr.h"

long
mpr_addmul_using_mpfr(mpr_ptr z, mpr_srcptr x, mpr_srcptr y, mp_bitcnt_t prec, mpfr_rnd_t rnd)
{
    long error;
    mpr_t t;

    mpr_init(t);
    mpr_mul_exact(t, x, y);
    error = mpr_add_using_mpfr(z, z, t, prec, rnd);
    mpr_clear(t);

    return error;
}

long
mpr_submul_using_mpfr(mpr_ptr z, mpr_srcptr x, mpr_srcptr y, mp_bitcnt_t prec, mpfr_rnd_t rnd)
{
    long error;
    mpr_t t;

    mpr_init(t);
    mpr_mul_exact(t, x, y);
    error = mpr_sub_using_mpfr(z, z, t, prec, rnd);
    mpr_clear(t);

    return error;
}

long
mpr_addmul_ui_using_mpfr(mpr_ptr z, mpr_srcptr x, ulong y, mp_bitcnt_t prec, mpfr_rnd_t rnd)
{
    long error;
    mpr_t t;
    mp_limb_t d;

    mpr_init(t);
    t->alloc = 1;
    t->d = &d;
    mpr_set_ui(t, y);

    error = mpr_addmul_using_mpfr(z, x, t, prec, rnd);

    return error;
}

long
mpr_submul_ui_using_mpfr(mpr_ptr z, mpr_srcptr x, ulong y, mp_bitcnt_t prec, mpfr_rnd_t rnd)
{
    long error;
    mpr_t t;
    mp_limb_t d;

    mpr_init(t);
    t->alloc = 1;
    t->d = &d;
    mpr_set_ui(t, y);

    error = mpr_submul_using_mpfr(z, x, t, prec, rnd);

    return error;
}
