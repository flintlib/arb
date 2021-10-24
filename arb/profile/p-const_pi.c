/*
    This is a significant rewrite of gmp-chudnovsky.c by Hanhong Xue with at
    least the following changes:

        bs_mul -> ui_factor_get_mpz_2exp uses squaring
        bs     -> pi_sum_split uses a basecase bigger than 1
        fmul   -> dust bin, there is no global factored number
        ...

    usage: make profile && ./build/arb/profile/p-const_pi

    known problems:
        - Memory usage is higher than expected.
        - Access to the global variable siever is not guarded by locks. Most
          accesses are read-only, but a write can only be prevented by fitting
          the table to the correct size at the beginning.
        - There is great confusion on how big mp_size_t is. It looks like it
          will overflow around a couple billion digits.
        - There is great confusion on what the return of mpn_mul is.
        - Untested with FLINT_BITS != 64
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <flint/profiler.h>
#include <flint/fmpz.h>
#include <flint/thread_pool.h>
#include <arb.h>

typedef unsigned int uint32_t;
#define MAXINT32 2147483647

/* TODO: try a proper implemention of unlikely */
#define unlikely(x) x
#define ASSERT(x) 
#define PROFILE_MEMORY 0

extern void arb_const_pi_chudnovsky_eval(arb_t s, slong prec);

/********************* factored integers *************************************/

typedef struct {
    ulong base, pow;
} ui_factor_entry;

// representation of an integer as prod_i base[i]^pow[i]

typedef struct {
    ui_factor_entry * data;
    slong alloc;
    ulong length;
} ui_factor_struct;

typedef ui_factor_struct ui_factor_t[1];


void ui_factor_init(ui_factor_t f)
{
    f->data = NULL;
    f->alloc = 0;
    f->length = 0;
}

void ui_factor_init2(ui_factor_t f, slong length)
{
    f->length = 0;
    f->alloc = length;
    if (length > 0)
        f->data = (ui_factor_entry *) flint_malloc(length * sizeof(ui_factor_entry));
    else
        f->data = NULL;
}

void ui_factor_clear(ui_factor_t f)
{
    if (f->data)
        flint_free(f->data);
}

/* the index of this struct represents the odd number n = 2*index + 1 */
/* index = 0 (n = 1) is special */
typedef struct {
  uint32_t pminus;      /* smallest prime factor p of n */
  uint32_t cofactor;    /* index (n/p - 1)/2 of the cofactor n/p */
} ui_factor_sieve_entry;

typedef struct {
    ui_factor_sieve_entry * array;
    ulong max;     /* the biggest number we can factor */
    slong alloc;
} ui_factor_sieve_struct;

typedef ui_factor_sieve_struct ui_factor_sieve_t[1];

static ui_factor_sieve_t siever;


void ui_factor_sieve_init(ui_factor_sieve_t S)
{
    S->array = NULL;
    S->alloc = 0;
    S->max = 0;
}

void ui_factor_sieve_clear(ui_factor_sieve_t S)
{
    #if PROFILE_MEMORY
        flint_printf("clearing ui_factor_sieve_t with alloc = %wd\n", S->alloc);
        flint_printf(" **** sieve used %wd KB\n", S->alloc*sizeof(ui_factor_sieve_entry)/1024);
    #endif
    if (S->array)
        flint_free(S->array);
}

void ui_factor_sieve_fit_length(ui_factor_sieve_t S, slong length)
{
    if (S->alloc >= length)
        return;

    S->array = (ui_factor_sieve_entry *) flint_realloc(S->array, length *
                                                sizeof(ui_factor_sieve_entry));
    S->alloc = length;
}

#define IDX(i) ((i)-1)/2

/* try to ensure that s can be used to factor numbers <= m */
void ui_factor_sieve_build(ui_factor_sieve_t S, ulong m)
{
    ulong i, j, k, n;

    m = FLINT_MIN(m, MAXINT32);
    m |= 1;

    if (m <= S->max)
        return;

    ui_factor_sieve_fit_length(S, (m + 1)/2);
    ui_factor_sieve_entry * s = S->array;

    s[IDX(1)].pminus = 0;
    s[IDX(1)].cofactor = 0;

    /*
        S->max = 0  =>  n = 1
        S->max = 1  =>  n = 1
        S->max = 2  =>  n = 1
        S->max = 3  =>  n = 3
        S->max = 4  =>  n = 3
    */
    n = FLINT_MAX(S->max, 1);
    n -= 1;
    n |= 1;

    for (i = n + 2; i <= m; i += 2)
    {
        s[IDX(i)].pminus = 0;
        s[IDX(i)].cofactor = 0;
    }

    for (i = 3; i <= n; i += 2)
    {
        ASSERT(s[IDX(i)].pminus >= 3);
        if (s[IDX(i)].cofactor != 0)
            continue;

        /* i is prime. scan odd k*i in (n, m] */
        for (k = m/i, k -= !(k&1); k*i > n; k -= 2)
        {
            if (s[IDX(k*i)].pminus != 0)
                continue;
            s[IDX(k*i)].pminus = i;
            s[IDX(k*i)].cofactor = IDX(k);
        }
    }

    for (i = n + 2; i <= m; i += 2)
    {
        if (s[IDX(i)].pminus != 0)
            continue;

        /* i is prime. scan odd k*i >= i^2 with k*i in (n, m] */
        s[IDX(i)].pminus = i;
        s[IDX(i)].cofactor = IDX(1);
        for (k = i; k <= m/i; k += 2)
        {
            if (s[IDX(k*i)].pminus != 0)
                continue;
            s[IDX(k*i)].pminus = i;
            s[IDX(k*i)].cofactor = IDX(k);
        }
    }

    for (i = 3; i <= m; i += 2)
    {
        ASSERT(s[IDX(i)].pminus >= 3);
        ASSERT(i == s[IDX(i)].pminus * (2*s[IDX(i)].cofactor + 1));
    }

    S->max = m;
}



void ui_factor_fit_length(ui_factor_t f, ulong length)
{
    if (f->alloc >= length)
        return;

    f->data = (ui_factor_entry *) flint_realloc(f->data, length *
                                                      sizeof(ui_factor_entry));
    f->alloc = length;
}

void ui_factor_swap(ui_factor_t f, ui_factor_t g)
{
    ui_factor_t t;
    *t = *f;
    *f = *g;
    *g = *t;
}

void ui_factor_print(const ui_factor_t f)
{
    int first = 1;
    const ui_factor_entry * fd = f->data;
    ulong fi;
    for (fi = 0; fi < f->length; fi++)
    {
        if (!first)
            flint_printf(" * ");
        flint_printf("%wu^%wu", fd[fi].base, fd[fi].pow);
        first = 0;
    }
    if (first)
        flint_printf("1");
}

int ui_factor_is_canonical(const ui_factor_t f)
{
    const ui_factor_entry * fd = f->data;
    ulong fi;
    for (fi = 0; fi < f->length; fi++)
    {
        if (fd[fi].pow == 0)
            return 0;

        if (fi > 0 && fd[fi-1].base >= fd[fi].base)
            return 0;
    }
    return 1;
}


void ui_factor_one(ui_factor_t f)
{
    f->length = 0;
}

void ui_factor_push_factor(ui_factor_t f, ulong base, ulong pow)
{
    ui_factor_fit_length(f, f->length + 1);
    f->data[f->length].base = base;
    f->data[f->length].pow = pow;
    f->length++;
}


void ui_factor_push_ui_no_sieve(ui_factor_t f, ulong b)
{
    ulong fi = f->length;
    n_factor_t fac;

    ui_factor_fit_length(f, fi + 16);

    ui_factor_entry * fd = f->data;

    n_factor_init(&fac);
    n_factor(&fac, b, 1);

    for (int i = 0; i < fac.num; i++)
    {
        fd[fi].base = fac.p[i];
        fd[fi].pow = fac.exp[i];
        fi++;
    }

    f->length = fi;
}

/* f *= base */
void ui_factor_push_ui(ui_factor_t f, ulong b, const ui_factor_sieve_t S)
{
    ulong fi = f->length;

    ui_factor_fit_length(f, fi + 16);

    ui_factor_entry * fd = f->data;

    if ((b % 2) == 0)
    {
        ulong e = 0;
        do {
            e += 1;
            b /= 2;
        } while ((b % 2) == 0);

        fd[fi].base = 2;
        fd[fi].pow = e;
        fi++;
    }

    if (unlikely(b > S->max))
    {
        f->length = fi;
        ui_factor_push_ui_no_sieve(f, b);
        return;
    }

    ui_factor_sieve_entry * s = S->array;

    b = b/2;
    if (b > 0)
    {
        fd[fi].base = s[b].pminus;
        b = s[b].cofactor;
        fd[fi].pow = 1;
        fi++;

        while (b > 0)
        {
            ulong p = s[b].pminus;
            b = s[b].cofactor;
            if (fd[fi-1].base == p)
            {
                fd[fi-1].pow++;
            }
            else
            {
                fd[fi].base = p;
                fd[fi].pow = 1;
                fi++;
            }
        }
    }

    f->length = fi;
}

/*
void mpn_sqr(mp_ptr a, mp_srcptr b, mp_size_t n)
{
    mpn_mul(a, b, n, b, n);
}
*/

/* x should have at least 3*len limbs allocated */
static mp_size_t ui_product_get_mpn(
    mp_limb_t * x,
    const ulong * v,
    ulong len,
    ulong stride)
{
    mp_size_t xlen, ylen, zlen;
    mp_limb_t out;

    ASSERT(0 < len);

    if (len <= 16)
    {
        xlen = 1;
        x[0] = v[0];
        for (ulong i = 1; i < len; i++)
        {
            out = mpn_mul_1(x, x, xlen, v[stride*i]);
            if (out != 0)
            {
                ASSERT(xlen < len);
                x[xlen] = out;
                xlen++;
            }
        }
    }
    else
    {
        ylen = ui_product_get_mpn(x + len, v, len/2 + (len % 2), 2*stride);
        zlen = ui_product_get_mpn(x + len + ylen, v + stride, len/2, 2*stride);
        out = (ylen >= zlen) ? mpn_mul(x, x + len, ylen, x + len + ylen, zlen)
                             : mpn_mul(x, x + len + ylen, zlen, x + len, ylen);
        xlen = ylen + zlen - (out == 0);
    }
    ASSERT(xlen <= len);
    ASSERT(xlen > 0);
    ASSERT(x[xlen - 1] != 0);
    return xlen;
}


static slong ui_product_one(ulong * terms)
{
    terms[0] = 1;
    return 0;
}

static slong ui_product_push(ulong * terms, slong top, ulong b)
{
    ulong hi, lo;
    umul_ppmm(hi, lo, terms[top], b);
    if (hi == 0)
    {
        terms[top] = lo;
    }
    else
    {
        terms[++top] = b;
    }
    return top;
}

static inline mp_limb_t * flint_mpz_fit_length(mpz_ptr z, mp_size_t n)
{
    if (n > z->_mp_alloc)
    {
        return (mp_ptr) _mpz_realloc(z, FLINT_MAX(z->_mp_alloc + z->_mp_alloc/2, n));
    }
    else
    {
        return z->_mp_d;
    }
}

/* power of 2 on the first base if the first base is 2 */
ulong _ui_factor_get_2exp(const ui_factor_t f)
{
    if (0 < f->length && f->data[0].base == 2)
        return f->data[0].pow;
    else
        return 0;
}

/* x = f forgetting about the first factor if the base is 2 */
void _ui_factor_get_mpz(
    mpz_t x,
    const ui_factor_t f,
    ui_factor_t t)  /* need temp space to clobber instead of f */
{
    ulong fi, ti, tj;
    ulong tn, fn = f->length;
    ui_factor_entry * td;
    const ui_factor_entry * fd = f->data;
    ulong b, p;
    mp_limb_t * xd;
    mp_size_t xn;
    mp_limb_t * zd;
    mp_size_t zn;
    mpz_t z;
    mpz_t y[FLINT_BITS + 1];
    mp_size_t l;
    slong i;
    mp_limb_t out;
    slong top;
    ulong * terms;
    TMP_INIT;

    fi = 0;

    if (fi < fn && fd[fi].base == 2)
    {
        fi++;
    }

    if (fi >= fn)
    {
        mpz_set_ui(x, 1);
        return;
    }

    TMP_START;

    td = TMP_ARRAY_ALLOC(fn, ui_factor_entry);
    terms = TMP_ARRAY_ALLOC(fn, ulong);

    ui_factor_fit_length(t, fn);
    td = t->data;

    top = ui_product_one(terms);
    tj = 0;
    for ( ; fi < fn; fi++)
    {
        b = fd[fi].base;
        p = fd[fi].pow;

        if (p % 2)
            top = ui_product_push(terms, top, b);

        p = p/2;
        if (p > 0)
        {
            td[tj].base = b;
            td[tj].pow = p;
            tj++;
        }
    }
    tn = tj;

    mpz_init(z);
    for (i = 0; i <= FLINT_BITS; i++)
        mpz_init(y[i]);

    i = 0;
    y[i]->_mp_size = ui_product_get_mpn(
                        flint_mpz_fit_length(y[i], 3*(top + 1)), terms, top + 1, 1);
    l = y[i]->_mp_size;

    while (tn > 0)
    {
        top = ui_product_one(terms);
        tj = 0;
        for (ti = 0; ti < tn; ti++)
        {
            b = td[ti].base;
            p = td[ti].pow;

            if (p % 2)
                top = ui_product_push(terms, top, b);

            p = p/2;
            if (p > 0)
            {
                td[tj].base = b;
                td[tj].pow = p;
                tj++;
            }
        }
        tn = tj;

        i++;
        y[i]->_mp_size = ui_product_get_mpn(
                        flint_mpz_fit_length(y[i], 3*(top + 1)), terms, top + 1, 1);
        l += y[i]->_mp_size << i;
    }

    if (i > 0)
    {
        xd = flint_mpz_fit_length(x, l);
        zd = flint_mpz_fit_length(z, l);

        zn = 2*y[i]->_mp_size;
        mpn_sqr(zd, y[i]->_mp_d, y[i]->_mp_size);
        zn -= (zd[zn - 1] == 0);

        while (1)
        {
            i--;
            xn = zn + y[i]->_mp_size;
            out = (zn >= y[i]->_mp_size) ?
                    mpn_mul(xd, zd, zn, y[i]->_mp_d, y[i]->_mp_size) :
                    mpn_mul(xd, y[i]->_mp_d, y[i]->_mp_size, zd, zn);
            xn -= (out == 0);

            ASSERT(xn <= l);

            if (i <= 0)
                break;

            mpn_sqr(zd, xd, xn);
            zn = 2*xn - (zd[2*xn - 1] == 0);
        }

        x->_mp_size = xn;
    }
    else
    {
        ASSERT(i == 0);
        mpz_set(x, y[i]);
    }

    mpz_clear(z);
    for (i = 0; i <= FLINT_BITS; i++)
        mpz_clear(y[i]);

    TMP_END;
}

ulong ui_factor_get_mpz_2exp(mpz_t x, const ui_factor_t f)
{
    ulong e;
    ui_factor_t t;
    ui_factor_init(t);
    _ui_factor_get_mpz(x, f, t);
    ui_factor_clear(t);
    return _ui_factor_get_2exp(f);
}

ulong ui_factor_get_fmpz_2exp(
    fmpz_t x,
    const ui_factor_t f)
{
    ulong e = ui_factor_get_mpz_2exp(_fmpz_promote(x), f);
    _fmpz_demote_val(x);
    return e;
}

int ui_factor_equal_fmpz(const ui_factor_t f, const fmpz_t x)
{
    fmpz_t y;
    int res;

    fmpz_init(y);
    ulong e = ui_factor_get_fmpz_2exp(y, f);
    fmpz_mul_2exp(y, y, e);
    res = fmpz_equal(y, x);
    fmpz_clear(y);
    return res;
}

int ui_factor_equal_mpz(const ui_factor_t f, const mpz_t x)
{
    fmpz_t y;
    fmpz_init(y);
    mpz_ptr Y = _fmpz_promote(y);
    ulong e = ui_factor_get_mpz_2exp(Y, f);
    mpz_mul_2exp(Y, Y, e);
    int res = mpz_cmp(Y, x) == 0;
    _fmpz_demote_val(y);
    fmpz_clear(y);
    return res;
}

// (f, g) = (f, g)/gcd(f, g)
void ui_factor_remove_gcd(ui_factor_t f, ui_factor_t g)
{
    ASSERT(ui_factor_is_canonical(f));
    ASSERT(ui_factor_is_canonical(g));

    ulong fi, gi, fj, gj;
    ulong fn = f->length;
    ulong gn = g->length;
    ui_factor_entry * fd = f->data;
    ui_factor_entry * gd = g->data;

    fi = gi = 0;
    fj = gj = 0;
    while (fj < fn && gj < gn)
    {
        if (fd[fj].base == gd[gj].base)
        {
            if (fd[fj].pow < gd[gj].pow)
            {
                ASSERT(gi <= gj);
                gd[gi].base = gd[gj].base;
                gd[gi].pow = gd[gj].pow - fd[fj].pow;
                gi++;
            }
            else if (fd[fj].pow > gd[gj].pow)
            {
                ASSERT(fi <= fj);
                fd[fi].base = fd[fj].base;
                fd[fi].pow = fd[fj].pow - gd[gj].pow;
                fi++;
            }
            fj++;
            gj++;
        }
        else if (fd[fj].base < gd[gj].base)
        {
            ASSERT(fi <= fj);
            fd[fi].base = fd[fj].base;
            fd[fi].pow = fd[fj].pow;
            fi++;
            fj++;
        }
        else
        {
            ASSERT(gi <= gj);
            gd[gi].base = gd[gj].base;
            gd[gi].pow = gd[gj].pow;
            gi++;
            gj++;
        }
    }

    while (fj < fn)
    {
        ASSERT(fi <= fj);
        fd[fi].base = fd[fj].base;
        fd[fi].pow = fd[fj].pow;
        fi++;
        fj++;
    }

    while (gj < gn)
    {
        ASSERT(gi <= gj);
        gd[gi].base = gd[gj].base;
        gd[gi].pow = gd[gj].pow;
        gi++;
        gj++;
    }

    f->length = fi;
    g->length = gi;

    ASSERT(ui_factor_is_canonical(f));
    ASSERT(ui_factor_is_canonical(g));
}


/* z = f*g */
static void ui_factor_mul(ui_factor_t z, const ui_factor_t f, const ui_factor_t g)
{
    ASSERT(ui_factor_is_canonical(f));
    ASSERT(ui_factor_is_canonical(g));

    ulong fn = f->length;
    ulong gn = g->length;
    ui_factor_fit_length(z, fn + gn);
    const ui_factor_entry * fd = f->data;
    const ui_factor_entry * gd = g->data;
    ui_factor_entry * zd = z->data;

    ulong fi = 0, gi = 0, zi = 0;
    while (fi < fn && gi < gn)
    {
        if (fd[fi].base == gd[gi].base)
        {
            zd[zi].base = fd[fi].base;
            zd[zi].pow = fd[fi].pow + gd[gi].pow;
            fi++;
            gi++;
        }
        else if (fd[fi].base < gd[gi].base)
        {
            zd[zi].base = fd[fi].base;
            zd[zi].pow = fd[fi].pow;
            fi++;
        }
        else
        {
            zd[zi].base = gd[gi].base;
            zd[zi].pow = gd[gi].pow;
            gi++;
        }

        zi++;
    }

    while (fi < fn)
    {
        zd[zi].base = fd[fi].base;
        zd[zi].pow = fd[fi].pow;
        fi++;
        zi++;
    }

    while (gi < gn)
    {
        zd[zi].base = gd[gi].base;
        zd[zi].pow = gd[gi].pow;
        gi++;
        zi++;
    }

    z->length = zi;

    ASSERT(ui_factor_is_canonical(z));
}


/* f = f^pow */
void ui_factor_pow_inplace(ui_factor_t f, ulong pow)
{
    ui_factor_entry * fd = f->data;
    ulong fn = f->length;
    for (ulong fi = 0; fi < fn; fi++)
        fd[fi].pow *= pow;
}

/* f *= g^p, like addmul on the exponents */
void ui_factor_mulpow_inplace(ui_factor_t f, const ui_factor_t g, ulong p)
{
    ui_factor_entry * fd = f->data;
    const ui_factor_entry * gd = g->data;
    ulong fn = f->length;
    ulong gn = g->length;
    ulong gi;
    ulong fi = 0;
    for (gi = 0; gi < gn; gi++)
    {
        while (1)
        {
            if (unlikely(fi >= fn))
            {
                ui_factor_fit_length(f, fi + gn - gi);
                fd = f->data;
                for ( ; gi < gn; gi++, fi++)
                {
                    fd[fi].base = gd[gi].base;
                    fd[fi].pow = p*gd[gi].pow;
                }
                f->length = fi;
                return;
            }

            if (fd[fi].base >= gd[gi].base)
                break;

            fi++;
        }

        if (unlikely(fd[fi].base > gd[gi].base))
        {
            ui_factor_fit_length(f, fn + 1);
            fd = f->data;
            ui_factor_entry t1, t2;
            t2.base = gd[gi].base;
            t2.pow = p*gd[gi].pow;
            for (ulong i = fi; i <= fn; i++)
            {
                t1 = fd[fi];
                fd[fi] = t2;
                t2 = t1;
            }
            fn++;
            fi++;
            f->length = fn;
        }
        else
        {
            fd[fi].pow += p*gd[gi].pow;
            fi++;
        }
    }
}

static void _ui_factor_sort_terms(ui_factor_entry * a,
                        slong length, ulong posmask, ulong totalmask)
{
    slong cur, mid;

    if (length < 20)
    {
        for (slong i = 1; i < length; i++)
        {
            for (slong j = i; j > 0 && a[j-1].base > a[j].base; j--)
            {
                ulong t1 = a[j].base;
                ulong t2 = a[j].pow;
                a[j].base =  a[j-1].base;
                a[j].pow = a[j-1].pow;
                a[j-1].base = t1;
                a[j-1].pow = t2;                
            }
        }
        return;
    }

    if ((posmask & totalmask) == 0)
    {
        posmask >>= 1;
        if (posmask != 0)
        {
            _ui_factor_sort_terms(a, length, posmask, totalmask);
        }

        return;
    }

    mid = 0;
    while (mid < length && (a[mid].base & posmask) == 0)
    {
        mid++;
    }

    cur = mid;
    while (++cur < length)
    {
        if ((a[cur].base & posmask) == 0)
        {
            ulong t1 = a[cur].base;
            ulong t2 = a[cur].pow;
            a[cur].base =  a[mid].base;
            a[cur].pow = a[mid].pow;
            a[mid].base = t1;
            a[mid].pow = t2;
            mid++;
        }
    }

    posmask >>= 1;
    if (posmask != 0)
    {
        _ui_factor_sort_terms(a, mid, posmask, totalmask);
        _ui_factor_sort_terms(a + mid, length - mid, posmask, totalmask);
    }
}

void ui_factor_canonicalise(ui_factor_t f)
{
    ui_factor_entry * fd = f->data;
    slong fn = f->length;
    slong i, j;

    if (unlikely(fn < 2))
        return;

    ulong smallest_idx = 0;
    ulong smallest = fd[smallest_idx].base;

    ulong base_mask = smallest;
    j = 1;
    for (i = 1; i < fn; i++)
    {
        ulong t1 = fd[i].base;
        ulong t2 = fd[i].pow;
        base_mask |= t1;
        if (t1 == smallest)
        {
            fd[smallest_idx].pow += t2;
        }
        else
        {
            if (t1 < smallest)
            {
                smallest = t1;
                smallest_idx = j;
            }
            fd[j].base = t1;
            fd[j].pow = t2;
            j++;
        }
    }

    fn = j;

    _ui_factor_sort_terms(fd, fn, UWORD(1) << (FLINT_BIT_COUNT(base_mask)-1), base_mask);

    j = 0;
    for (i = 1; i < fn; i++)
    {
        ASSERT(j < i);

        if (fd[j].base == fd[i].base)
        {
            fd[j].pow += fd[i].pow;
        }
        else
        {
            j++;
            fd[j].base = fd[i].base;
            fd[j].pow = fd[i].pow;
        }
    }

    j++;
    f->length = j;
}

mp_size_t flint_mpn_mul_11(mp_limb_t * y, const mp_limb_t * x, mp_size_t n, mp_limb_t a1, mp_limb_t a2)
{
    mp_limb_t hi, lo;

    umul_ppmm(hi, lo, a1, a2);
    if (hi == 0)
    {
        y[n] = mpn_mul_1(y, x, n, lo); n += (y[n] != 0);
        return n;
    }

    y[n] = mpn_mul_1(y, x, n, a1); n += (y[n] != 0);
    y[n] = mpn_mul_1(y, y, n, a2); n += (y[n] != 0);
    return n;
}


mp_size_t flint_mpn_mul_111(mp_limb_t * y, const mp_limb_t * x, mp_size_t n, mp_limb_t a1, mp_limb_t a2, mp_limb_t a3)
{
    mp_limb_t hi, lo, p1, p2;

    umul_ppmm(hi, lo, a1, a2);
    if (hi == 0)
    {
        p1 = lo;
        umul_ppmm(hi, lo, p1, a3);
        if (hi == 0)
        {
            y[n] = mpn_mul_1(y, x, n, lo); n += (y[n] != 0);
            return n;
        }
        p2 = a3;
        goto do_two;
    }

    umul_ppmm(hi, lo, a1, a3);
    if (hi == 0)
    {
        p1 = lo;
        p2 = a2;
        goto do_two;
    }

    umul_ppmm(hi, lo, a2, a3);
    if (hi == 0)
    {
        p1 = lo;
        p2 = a1;
        goto do_two;
    }    

    y[n] = mpn_mul_1(y, x, n, a1); n += (y[n] != 0);
    y[n] = mpn_mul_1(y, y, n, a2); n += (y[n] != 0);
    y[n] = mpn_mul_1(y, y, n, a3); n += (y[n] != 0);
    return n;

do_two:

    y[n] = mpn_mul_1(y, x, n, p1); n += (y[n] != 0);
    y[n] = mpn_mul_1(y, y, n, p2); n += (y[n] != 0);
    return n;
}

/**************************** start pi stuff *********************************/
/* my result (ON CURRENT BOX!!!):
virt/peak/res/peak(MB): 21.12 336.54 6.35 272.19
****** one thread ********
prec      1024: here      1, arb      1, ratio 1.00    diff*2^prec: 0.35767 +/- 3.1035
prec      2048: here      1, arb      1, ratio 1.00    diff*2^prec: -0.20294 +/- 3.6166
prec      4096: here      1, arb      1, ratio 1.00    diff*2^prec: 0.8522 +/- 3.1691
prec      8192: here      1, arb      1, ratio 1.00    diff*2^prec: 0.39326 +/- 2.4516
prec     16384: here      1, arb      1, ratio 1.00    diff*2^prec: 0.14836 +/- 2.3905
prec     32768: here      1, arb      1, ratio 1.00    diff*2^prec: 0.53247 +/- 2.1235
prec     65536: here      1, arb      2, ratio 2.00    diff*2^prec: 0.19283 +/- 2.234
prec    131072: here      2, arb      4, ratio 2.00    diff*2^prec: 0.88084 +/- 2.1784
prec    262144: here      7, arb      9, ratio 1.29    diff*2^prec: 0.31244 +/- 2.3353
prec    524288: here     16, arb     24, ratio 1.50    diff*2^prec: 0.15431 +/- 2.3108
prec   1048576: here     38, arb     55, ratio 1.45    diff*2^prec: 0.059249 +/- 2.2084
prec   2097152: here     89, arb    133, ratio 1.49    diff*2^prec: 0.14183 +/- 2.3491
prec   4194304: here    206, arb    314, ratio 1.52    diff*2^prec: -0.054642 +/- 2.1233
prec   8388608: here    485, arb    741, ratio 1.53    diff*2^prec: 0.64171 +/- 2.0171
prec  16777216: here   1149, arb   1760, ratio 1.53    diff*2^prec: 0.7042 +/- 2.1601
prec  33554432: here   2727, arb   4192, ratio 1.54    diff*2^prec: 1.0175 +/- 2.13
prec  67108864: here   6378, arb  10006, ratio 1.57    diff*2^prec: 0.64759 +/- 2.0933
prec 134217728: here  14756, arb  23442, ratio 1.59    diff*2^prec: 0.21221 +/- 2.1246
****** two threads *******
prec      1024: here      1, arb      1, ratio 1.00    diff*2^prec: 0.85767 +/- 2.4476
prec      2048: here      1, arb      1, ratio 1.00    diff*2^prec: 0.29706 +/- 2.3662
prec      4096: here      1, arb      1, ratio 1.00    diff*2^prec: 0.8522 +/- 2.1024
prec      8192: here      1, arb      1, ratio 1.00    diff*2^prec: 0.39326 +/- 2.2593
prec     16384: here      1, arb      1, ratio 1.00    diff*2^prec: 0.64836 +/- 2.2141
prec     32768: here      1, arb      1, ratio 1.00    diff*2^prec: 0.032475 +/- 2.0377
prec     65536: here      1, arb      2, ratio 2.00    diff*2^prec: 0.19283 +/- 2.2527
prec    131072: here      2, arb      4, ratio 2.00    diff*2^prec: 0.38084 +/- 2.1352
prec    262144: here      5, arb      9, ratio 1.80    diff*2^prec: 0.31244 +/- 2.3356
prec    524288: here     11, arb     24, ratio 2.18    diff*2^prec: 0.15431 +/- 2.1656
prec   1048576: here     26, arb     57, ratio 2.19    diff*2^prec: 0.55925 +/- 2.461
prec   2097152: here     61, arb    135, ratio 2.21    diff*2^prec: 0.14183 +/- 2.3491
prec   4194304: here    140, arb    319, ratio 2.28    diff*2^prec: 0.94536 +/- 2.278
prec   8388608: here    326, arb    748, ratio 2.29    diff*2^prec: 1.1417 +/- 2.279
prec  16777216: here    760, arb   1770, ratio 2.33    diff*2^prec: 0.7042 +/- 2.1601
prec  33554432: here   1801, arb   4229, ratio 2.35    diff*2^prec: 0.51755 +/- 2.3609
prec  67108864: here   4146, arb  10081, ratio 2.43    diff*2^prec: 0.64759 +/- 2.0933
prec 134217728: here   9538, arb  23622, ratio 2.48    diff*2^prec: 0.71221 +/- 2.3141
*/

/* timing break down in ms for the 10^6 digit computation (ON AN OLD BOX!!!):
    485ms total
     = 53ms basecase
     + 269ms fold
       =  85ms ui_factor_get_fmpz_2exp
       +  12ms ui_facor mul/gcd
       + 171ms fmpz add/mul
     + 150ms float at the end
*/


/*
    Set (p, r, q) = sum of term_j with j in [start,stop] inclusive.
    The factored number mult is 640320^3/24.
    s is temporary working space
*/
static void pi_sum_basecase(
    mpz_t p, mpz_t s, ui_factor_t r, ui_factor_t q,
    ulong start, ulong stop,
    const ui_factor_t mult)
{
    mp_limb_t out, a[2];
    mp_size_t p_n, s_n;
    mp_limb_t * p_d, * s_d;
    ulong j;

    ASSERT(1 <= start);
    ASSERT(start <= stop);

    /* good guesses on how big things will get */
    ui_factor_fit_length(r, 8*(stop - start + 1));
    ui_factor_fit_length(q, 4*(stop - start + 1));
    p_d = flint_mpz_fit_length(p, 2*(stop - start) + 20);
    s_d = flint_mpz_fit_length(s, 1*(stop - start) + 10);

    ui_factor_one(r);
    ui_factor_one(q);

    j = start;

    // q = j;
    ui_factor_push_ui(q, j, siever);

    // r = s = (2*j-1)*(6*j-1)*(6*j-5), s will be expanded form of r
    ui_factor_push_ui(r, 2*j-1, siever);
    ui_factor_push_ui(r, 6*j-1, siever);
    ui_factor_push_ui(r, 6*j-5, siever);

    s_d[0] = 6*j-1;
    s_n = flint_mpn_mul_11(s_d, s_d, 1, 6*j-5, 2*j-1);

    // a = 13591409 + 545140134*j
    umul_ppmm(a[1], a[0], j, 545140134);
    add_ssaaaa(a[1], a[0], a[1], a[0], 0, 13591409);

    // p = s*a
    if (a[1] == 0)
    {
        p_d[s_n] = mpn_mul_1(p_d, s_d, s_n, a[0]);
        p_n = s_n + (p_d[s_n] != 0);
    }
    else
    {
        out = mpn_mul_1(p_d + 1, s_d, s_n, a[1]);
        p_d[s_n + 1] = out;
        p_n = s_n + 1 + (out != 0);

        out = mpn_addmul_1(p_d, s_d, s_n, a[0]);
        if (out != 0)
            out = mpn_add_1(p_d + s_n, p_d + s_n, p_n - s_n, out);
        p_d[p_n] = out;
        p_n += (out != 0);
    }

    for (++j; j <= stop; ++j)
    {
        // q *= j
        ui_factor_push_ui(q, j, siever);

        // r, s *= (2*j-1)*(6*j-1)*(6*j-5)
        ui_factor_push_ui(r, 2*j-1, siever);
        ui_factor_push_ui(r, 6*j-5, siever);
        ui_factor_push_ui(r, 6*j-1, siever);
        s_n = flint_mpn_mul_111(s_d, s_d, s_n, 2*j-1, 6*j-5, 6*j-1);

        // p *= 10939058860032000*j^3
        p_n = flint_mpn_mul_111(p_d, p_d, p_n, j, j, j);
        if (FLINT_BITS == 64)
        {
            p_d[p_n] = mpn_mul_1(p_d, p_d, p_n, UWORD(10939058860032000)); p_n += (p_d[p_n] != 0);
        }
        else
        {
            p_d[p_n] = mpn_mul_1(p_d, p_d, p_n, 296740963); p_n += (p_d[p_n] != 0);
            p_d[p_n] = mpn_mul_1(p_d, p_d, p_n, 36864000); p_n += (p_d[p_n] != 0);
        }

        // a = 13591409 + 545140134*j
        add_ssaaaa(a[1], a[0], a[1], a[0], 0, 545140134);

        // p +-= s*a
        ASSERT(p_n >= s_n);
        if ((j & 1) ^ (start & 1))
        {
            if (a[1] != 0)
            {
                out = mpn_submul_1(p_d + 1, s_d, s_n, a[1]);
                if (p_n > s_n + 1 && out != 0)
                    out = mpn_sub_1(p_d + (s_n + 1), p_d + (s_n + 1), p_n - (s_n + 1), out);
                ASSERT(out == 0);
                while (p_d[p_n - 1] == 0)
                {
                    ASSERT(p_n > 0);
                    p_n--;
                }
            }

            ASSERT(p_n >= s_n);

            out = mpn_submul_1(p_d, s_d, s_n, a[0]);
            if (p_n > s_n && out != 0)
                out = mpn_sub_1(p_d + s_n, p_d + s_n, p_n - s_n, out);
            ASSERT(out == 0);
            while (p_d[p_n - 1] == 0)
            {
                ASSERT(p_n > 0);
                p_n--;
            }
        }
        else
        {
            out = mpn_addmul_1(p_d, s_d, s_n, a[0]);
            if (p_n > s_n && out != 0)
                out = mpn_add_1(p_d + s_n, p_d + s_n, p_n - s_n, out);
            p_d[p_n] = out;
            p_n += (out != 0);

            if (a[1] != 0)
            {
                ASSERT(p_n >= s_n + 1);                
                out = mpn_addmul_1(p_d + 1, s_d, s_n, a[1]);
                if (p_n > s_n + 1 && out != 0)
                    out = mpn_add_1(p_d + s_n + 1, p_d + s_n + 1, p_n - (s_n + 1), out);
                p_d[p_n] = out;
                p_n += (out != 0);
            }
        }

        s_d = flint_mpz_fit_length(s, s_n + 3);
        p_d = flint_mpz_fit_length(p, p_n + 4);
    }

    s->_mp_size = s_n;
    p->_mp_size = (start & 1) ? -p_n : p_n;

    ui_factor_canonicalise(q);
    ui_factor_pow_inplace(q, 3);
    ui_factor_mulpow_inplace(q, mult, stop - start + 1);

    ui_factor_canonicalise(r);

    ASSERT(ui_factor_equal_mpz(r, s));
}

// (p1, r1, q1) = (p1*q2 + r1*p2, r1*r2, q1*q2) / gcd(r1, q2)
static void fold(
    mpz_t p1, ui_factor_t r1, ui_factor_t q1,
    mpz_t p2, ui_factor_t r2, ui_factor_t q2,
    int needr,  /* bool */
    mpz_t t2,   /* temps */
    mpz_t t1,
    ui_factor_t s)
{
    ulong q2e, r1e;

    ui_factor_remove_gcd(r1, q2);

    // p1 = p1*q2 + r1*p2

    _ui_factor_get_mpz(t2, q2, s);
    mpz_mul(t1, p1, t2);
    _ui_factor_get_mpz(t2, r1, s);
    mpz_mul(p1, p2, t2);
    q2e = _ui_factor_get_2exp(q2);
    r1e = _ui_factor_get_2exp(r1);
    if (q2e >= r1e)
        mpz_mul_2exp(t1, t1, q2e - r1e);
    else
        mpz_mul_2exp(p1, p1, r1e - q2e);
    mpz_add(p1, p1, t1);
    mpz_mul_2exp(p1, p1, FLINT_MIN(r1e, q2e));

    // q1 = q1*q2
    ui_factor_mul(s, q1, q2);
    ui_factor_swap(q1, s);

    if (needr)
    {
        // r1 = r1*r2
        ui_factor_mul(s, r1, r2);
        ui_factor_swap(r1, s);
    }
}

// Set (p, r, q) = sum of terms in [start,stop] inclusive.
static void pi_sum_split(
    mpz_t p, ui_factor_t r, ui_factor_t q,
    ulong start, ulong stop,
    int needr, /* bool */
    const ui_factor_t mult)
{
    mpz_t p1, t1, t2;
    ui_factor_t r1, q1, s;
    ulong diff, mid;

    ASSERT(start <= stop);

    mpz_init(p1);
    ui_factor_init(r1);
    ui_factor_init(q1);

    mpz_init(t1);
    mpz_init(t2);
    ui_factor_init(s);

    diff = stop - start;
    if (diff > 300)
    {
        mid = diff/16*9 + start;
        pi_sum_split(p, r, q, start, mid, 1, mult);
        pi_sum_split(p1, r1, q1, mid + 1, stop, 1, mult);
        fold(p, r, q, p1, r1, q1, needr, t1, t2, s);
    }
    else if (diff > 200)
    {
        ulong mid1 = diff/4 + start;
        ulong mid2 = diff/2 + start;
        ulong mid3 = diff/4 + diff/2 + start;

        pi_sum_basecase(p, t1, r, q, start, mid1, mult);
        pi_sum_basecase(p1, t1, r1, q1, mid1+1, mid2, mult);
        fold(p, r, q, p1, r1, q1, 1, t1, t2, s);
        pi_sum_basecase(p1, t1, r1, q1, mid2+1, mid3, mult);
        fold(p, r, q, p1, r1, q1, 1, t1, t2, s);
        pi_sum_basecase(p1, t1, r1, q1, mid3+1, stop, mult);
        fold(p, r, q, p1, r1, q1, 1, t1, t2, s);
    }
    else if (diff > 100)
    {
        mid = diff/2 + start;
        pi_sum_basecase(p, t1, r, q, start, mid, mult);
        pi_sum_basecase(p1, t2, r1, q1, mid + 1, stop, mult);
        fold(p, r, q, p1, r1, q1, needr, t1, t2, s);
    }
    else
    {
        pi_sum_basecase(p, t1, r, q, start, stop, mult);
    }

    mpz_clear(p1);
    ui_factor_clear(r1);
    ui_factor_clear(q1);

    mpz_clear(t1);
    mpz_clear(t2);
    ui_factor_clear(s);
}


typedef struct {
    mpz_t p2;
    ui_factor_t r2;
    ui_factor_t q2;
    ui_factor_struct * mult;
    ulong start;
    ulong stop;
    ulong q2e;
    mpz_t p1q2;
    mpz_ptr p1;
} worker_arg;

void worker_proc(void * varg)
{
    worker_arg * arg = (worker_arg *) varg;
    pi_sum_split(arg->p2, arg->r2, arg->q2, arg->start, arg->stop, 1, arg->mult);
}

void worker_proc2(void * varg)
{
    mpz_t q2t;
    worker_arg * arg = (worker_arg *) varg;

    mpz_init(q2t);
    arg->q2e = ui_factor_get_mpz_2exp(q2t, arg->q2);
    mpz_mul(arg->p1q2, arg->p1, q2t);
    mpz_clear(q2t);
}


ulong pi_sum(fmpz_t p, fmpz_t q, ulong num_terms)
{
    ulong qe;
    ui_factor_t r1, q1, mult;
    mpz_ptr p1 = _fmpz_promote(p);

    ui_factor_init(r1);
    ui_factor_init(q1);
    ui_factor_init(mult);

    ui_factor_sieve_init(siever);
    ui_factor_sieve_build(siever, FLINT_MAX(UWORD(3*5*23*29), 6*num_terms-1));

    ui_factor_one(mult);
    ui_factor_push_factor(mult, 2, 15);
    ui_factor_push_factor(mult, 3, 2);
    ui_factor_push_factor(mult, 5, 3);
    ui_factor_push_factor(mult, 23, 3);
    ui_factor_push_factor(mult, 29, 3);

    if (global_thread_pool_initialized && num_terms > 20)
    {
        thread_pool_handle handles[1];
        slong num_handles = thread_pool_request(global_thread_pool, handles, 1);
        if (num_handles > 0)
        {
            /*
                This can probably be cleaned up alot. Since the floating point
                arithmetic was not threaded at the time of writing, there was/is
                little point in going past two threads.
            */
            ulong mid = num_terms/16*9;
            ulong r1e;
            mpz_t r1t, p2r1;
            ui_factor_t q1q2;
            worker_arg warg;

            ASSERT(num_handles == 1);

            mpz_init(warg.p2);
            mpz_init(warg.p1q2);
            ui_factor_init(warg.r2);
            ui_factor_init(warg.q2);

            warg.mult = mult;
            warg.p1 = p1;
            warg.start = mid + 1;
            warg.stop = num_terms;

            /* calculate [1, mid] and [mid + 1, num_terms] */

            thread_pool_wake(global_thread_pool, handles[0], 1, &worker_proc, &warg);
            pi_sum_split(p1, r1, q1, 1, mid, 1, mult);
            thread_pool_wait(global_thread_pool, handles[0]);

            /* join the two pieces */

            ui_factor_remove_gcd(r1, warg.q2);

            thread_pool_wake(global_thread_pool, handles[0], 1, &worker_proc2, &warg);

            mpz_init(r1t);
            r1e = ui_factor_get_mpz_2exp(r1t, r1);

            mpz_init(p2r1);
            mpz_mul(p2r1, warg.p2, r1t);

            ui_factor_init(q1q2);
            ui_factor_mul(q1q2, q1, warg.q2);
            qe = ui_factor_get_fmpz_2exp(q, q1q2);

            thread_pool_wait(global_thread_pool, handles[0]);

            if (warg.q2e >= r1e)
            {
                mpz_mul_2exp(p1, warg.p1q2, warg.q2e - r1e);
                mpz_add(p1, p1, p2r1);
                mpz_mul_2exp(p1, p1, r1e);
            }
            else
            {
                mpz_mul_2exp(p1, p2r1, r1e - warg.q2e);
                mpz_add(p1, p1, warg.p1q2);
                mpz_mul_2exp(p1, p1, warg.q2e);
            }

            thread_pool_give_back(global_thread_pool, handles[0]);

            ui_factor_clear(q1q2);
            mpz_clear(r1t);
            mpz_clear(p2r1);

            mpz_clear(warg.p2);
            mpz_clear(warg.p1q2);
            ui_factor_clear(warg.r2);
            ui_factor_clear(warg.q2);

            goto cleanup;
        }
    }

    /* serial */
    pi_sum_split(p1, r1, q1, 1, num_terms, 0, mult);
    qe = ui_factor_get_fmpz_2exp(q, q1);

cleanup:

    ui_factor_sieve_clear(siever);

    ui_factor_clear(r1);
    ui_factor_clear(q1);
    ui_factor_clear(mult);

    _fmpz_demote_val(p);

    return qe;
}

void pi_here(arb_t U, slong prec)
{
    slong wp = prec + 3;
    ulong num_terms = prec * 0.021226729578153557 + 2;
    ulong qe;
    fmpz_t p, q;
    arb_t P, Q, T;

    fmpz_init(p);
    fmpz_init(q);
    arb_init(P);
    arb_init(Q);
    arb_init(T);

    qe = pi_sum(p, q, num_terms);
    /*
        we now have p/(q*2^qe) = sum from 1 to num_terms
        finish off the pi calculation with

                 q*2^qe * 640320/12
        -----------------------------------
        (p + 13591409*q*2^qe)*rsqrt(640320)
    */
    arb_set_round_fmpz(Q, q, wp);
    arb_mul_2exp_si(Q, Q, qe);
    arb_mul_ui(Q, Q, 640320/12, wp);
    fmpz_mul_ui(q, q, 13591409);
    fmpz_mul_2exp(q, q, qe);
    fmpz_add(q, q, p);
    arb_set_round_fmpz(P, q, wp);
    arb_rsqrt_ui(U, 640320, wp);
    arb_mul(T, P, U, wp);
    arb_div(U, Q, T, wp);

    fmpz_clear(p);
    fmpz_clear(q);
    arb_clear(P);
    arb_clear(Q);
    arb_clear(T);
}

void compare_pi(slong prec)
{
    timeit_t timer;
    slong t1, t2;
    arb_t T, U, U2;

    arb_init(T);
    arb_init(U);
    arb_init(U2);

    flint_printf("prec %9wd: ", prec);
    fflush(stdout);

    timeit_start(timer);
    pi_here(U, prec);
    timeit_stop(timer);
    t1 = timer->wall;
    t1 = FLINT_MAX(WORD(1), t1);

    flint_printf("here %6wd, ", t1);
    fflush(stdout);

    timeit_start(timer);
    arb_const_pi_chudnovsky_eval(U2, prec);
    timeit_stop(timer);
    t2 = timer->wall;

    t2 = FLINT_MAX(WORD(1), t2);

    flint_printf("arb %6wd, ratio %0.2f", t2, (double)t2/(double)t1);

    arb_sub(T, U, U2, prec);
    arb_mul_2exp_si(T, T, prec);
    flint_printf("    diff*2^prec: ");
    arb_printd(T, 5);
    flint_printf("\n");
    fflush(stdout);

    arb_clear(T);
    arb_clear(U);
    arb_clear(U2);
}


int main(int i, char * b)
{
    if (1)
    {
        arb_t u;
        arb_init(u);
        pi_here(u, WORD(100000000));
        arb_clear(u);
        SHOW_MEMORY_USAGE;
    }

    printf("****** one thread ********\n");
    for (i = 10; i < 28; i++)
        compare_pi(WORD(1)<<i);

    printf("****** two threads *******\n");
    flint_set_num_threads(2);
    for (i = 10; i < 28; i++)
        compare_pi(WORD(1)<<i);

    flint_cleanup_master();
}

