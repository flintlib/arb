/*
    Copyright (C) 2017 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "mag.h"

unsigned int randman(flint_rand_t state)
{
    switch (n_randint(state, 8))
    {
        case 0:
        case 1:
        case 2:
        case 3:
            return (1 << (MAG_BITS)) - 1;
        case 4:
        case 5:
            return (1 << (MAG_BITS - 1));
        default:
            return (n_randlimb(state) % (1 << MAG_BITS)) | (1 << (MAG_BITS - 1));
    }
}

slong randexp(flint_rand_t state)
{
    switch (n_randint(state, 8))
    {
        case 0:
        case 1:
        case 2:
        case 3:
            return 0;
        case 4:
        case 5:
            return 1;
        default:
            return n_randint(state, 2 * MAG_BITS) - MAG_BITS;
    }
}

void test_MAG_ADD()
{
    slong iter;
    flint_rand_t state;

    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000000 * arb_test_multiplier(); iter++)
    {
        unsigned int a, b, r;
        double aa, bb, rr;
        slong aexp, bexp, rexp;

        a = randman(state);
        b = randman(state);
        r = randman(state);

        aexp = randexp(state);
        bexp = randexp(state);
        rexp = randexp(state);

        aa = ldexp(a, aexp - MAG_BITS);
        bb = ldexp(b, bexp - MAG_BITS);

        MAG_ADD(r, rexp, a, aexp, b, bexp);

        rr = ldexp(r, rexp - MAG_BITS);

        if (FLINT_BIT_COUNT(r) != MAG_BITS)
        {
            flint_printf("fail: mantissa has %d bits!\n", FLINT_BIT_COUNT(r));
            flint_abort();
        }

        if (rr < (aa + bb) * (1 - 1e-14))
        {
            flint_printf("fail: result is too small!\n");
            flint_abort();
        }

        if (rr > (aa + bb) * (1 + 1e-6))
        {
            flint_printf("fail: result is too large!\n");
            flint_printf("aa = %f\n", aa);
            flint_printf("bb = %f\n", bb);
            flint_printf("rr = %f\n", rr);
            flint_printf("aa + bb = %f\n", aa + bb);
            flint_abort();
        }
    }

    flint_randclear(state);
}

void test_MAG_ADD_2EXP()
{
    slong iter;
    flint_rand_t state;

    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000000 * arb_test_multiplier(); iter++)
    {
        unsigned int a, r;
        double aa, bb, rr;
        slong aexp, bexp, rexp;

        a = randman(state);
        r = randman(state);

        aexp = randexp(state);
        bexp = randexp(state);
        rexp = randexp(state);

        aa = ldexp(a, aexp - MAG_BITS);
        bb = ldexp(1.0, bexp);

        MAG_ADD_2EXP(r, rexp, a, aexp, bexp);

        rr = ldexp(r, rexp - MAG_BITS);

        if (FLINT_BIT_COUNT(r) != MAG_BITS)
        {
            flint_printf("fail: mantissa has %d bits!\n", FLINT_BIT_COUNT(r));
            flint_abort();
        }

        if (rr < (aa + bb) * (1 - 1e-14))
        {
            flint_printf("fail: result is too small!\n");
            flint_abort();
        }

        if (rr > (aa + bb) * (1 + 1e-6))
        {
            flint_printf("fail: result is too large!\n");
            flint_printf("aa = %f\n", aa);
            flint_printf("bb = %f\n", bb);
            flint_printf("rr = %f\n", rr);
            flint_printf("aa + bb = %f\n", aa + bb);
            flint_abort();
        }
    }

    flint_randclear(state);
}

void test_MAG_MUL()
{
    slong iter;
    flint_rand_t state;

    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000000 * arb_test_multiplier(); iter++)
    {
        unsigned int a, b, r;
        double aa, bb, rr;
        slong aexp, bexp, rexp;

        a = randman(state);
        b = randman(state);
        r = randman(state);

        aexp = randexp(state);
        bexp = randexp(state);
        rexp = randexp(state);

        aa = ldexp(a, aexp - MAG_BITS);
        bb = ldexp(b, bexp - MAG_BITS);

        MAG_MUL(r, rexp, a, aexp, b, bexp);

        rr = ldexp(r, rexp - MAG_BITS);

        if (FLINT_BIT_COUNT(r) != MAG_BITS)
        {
            flint_printf("fail: mantissa has %d bits!\n", FLINT_BIT_COUNT(r));
            flint_abort();
        }

        if (rr < (aa * bb) * (1 - 1e-14))
        {
            flint_printf("fail: result is too small!\n");
            flint_abort();
        }

        if (rr > (aa * bb) * (1 + 1e-6))
        {
            flint_printf("fail: result is too large!\n");
            flint_printf("aa = %f\n", aa);
            flint_printf("bb = %f\n", bb);
            flint_printf("rr = %f\n", rr);
            flint_printf("aa * bb = %f\n", aa * bb);
            flint_abort();
        }
    }

    flint_randclear(state);
}

void test_MAG_ADDMUL()
{
    slong iter;
    flint_rand_t state;

    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000000 * arb_test_multiplier(); iter++)
    {
        unsigned int a, b, c, r;
        double aa, bb, cc, rr;
        slong aexp, bexp, cexp, rexp;

        a = randman(state);
        b = randman(state);
        c = randman(state);
        r = randman(state);

        aexp = randexp(state);
        bexp = randexp(state);
        cexp = randexp(state);
        rexp = randexp(state);

        aa = ldexp(a, aexp - MAG_BITS);
        bb = ldexp(b, bexp - MAG_BITS);
        cc = ldexp(c, cexp - MAG_BITS);

        MAG_ADDMUL(r, rexp, a, aexp, b, bexp, c, cexp);

        rr = ldexp(r, rexp - MAG_BITS);

        if (FLINT_BIT_COUNT(r) != MAG_BITS)
        {
            flint_printf("fail: mantissa has %d bits!\n", FLINT_BIT_COUNT(r));
            flint_abort();
        }

        if (rr < (aa + bb * cc) * (1 - 1e-14))
        {
            flint_printf("fail: result is too small!\n");
            flint_abort();
        }

        if (rr > (aa + bb * cc) * (1 + 1e-6))
        {
            flint_printf("fail: result is too large!\n");
            flint_printf("aa = %f\n", aa);
            flint_printf("bb = %f\n", bb);
            flint_printf("cc = %f\n", cc);
            flint_printf("rr = %f\n", rr);
            flint_printf("aa + bb * cc = %f\n", aa + bb * cc);
            flint_abort();
        }
    }

    flint_randclear(state);
}

void
test_MAG_MULADDMUL()
{
    slong iter;
    flint_rand_t state;

    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000000 * arb_test_multiplier(); iter++)
    {
        unsigned int a, b, c, d, r;
        double aa, bb, cc, dd, rr;
        slong aexp, bexp, cexp, dexp, rexp;

        a = randman(state);
        b = randman(state);
        c = randman(state);
        d = randman(state);
        r = randman(state);

        aexp = randexp(state);
        bexp = randexp(state);
        cexp = randexp(state);
        dexp = randexp(state);
        rexp = randexp(state);

        aa = ldexp(a, aexp - MAG_BITS);
        bb = ldexp(b, bexp - MAG_BITS);
        cc = ldexp(c, cexp - MAG_BITS);
        dd = ldexp(d, dexp - MAG_BITS);

        MAG_MULADDMUL(r, rexp, a, aexp, b, bexp, c, cexp, d, dexp);

        rr = ldexp(r, rexp - MAG_BITS);

        if (FLINT_BIT_COUNT(r) != MAG_BITS)
        {
            flint_printf("fail: mantissa has %d bits!\n", FLINT_BIT_COUNT(r));
            flint_abort();
        }

        if (rr < (aa * bb + cc * dd) * (1 - 1e-14))
        {
            flint_printf("fail: result is too small!\n");
            flint_abort();
        }

        if (rr > (aa * bb + cc * dd) * (1 + 1e-6))
        {
            flint_printf("fail: result is too large!\n");
            flint_printf("aa = %f\n", aa);
            flint_printf("bb = %f\n", bb);
            flint_printf("cc = %f\n", cc);
            flint_printf("dd = %f\n", dd);
            flint_printf("rr = %f\n", rr);
            flint_printf("aa * bb + cc * dd = %f\n", aa * bb + cc * dd);
            flint_abort();
        }
    }

    flint_randclear(state);
}

int main()
{
    flint_printf("macros....");
    fflush(stdout);

    test_MAG_ADD();
    test_MAG_ADD_2EXP();
    test_MAG_MUL();
    test_MAG_ADDMUL();
    test_MAG_MULADDMUL();

    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

