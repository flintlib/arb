/*
    Copyright (C) 2016 Pascal Molin

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_dirichlet.h"

void
test_dft(ulong q)
{
    ulong i;
    slong prec = 100;
    acb_dirichlet_group_t G;
    acb_dirichlet_conrey_t x;
    acb_dirichlet_char_t chi;
    acb_t s, z;
    acb_ptr v;

    acb_dirichlet_group_init(G, q);
    acb_dirichlet_conrey_init(x, G);
    acb_dirichlet_char_init(chi, G);

    acb_init(s);
    acb_one(s);
    acb_div_si(s, s, 2, prec);

    v = _acb_vec_init(G->phi_q);

    /* all at once */
    acb_dirichlet_l_vec_hurwitz(v, s, G, prec);

    /* check with complete loop */

    i = 0;
    acb_init(z);
    acb_dirichlet_conrey_one(x, G);
    do {

        acb_dirichlet_char_conrey(chi, G, x);
        acb_dirichlet_l_hurwitz(z, s, G, chi, prec);

        if (!acb_overlaps(z, v + i))
        {
            flint_printf("\n L value differ");
            flint_printf("\nL(1/2, %wu) single = ", x->n);
            acb_printd(z, 20);
            flint_printf("\nL(1/2, %wu) multi = ", x->n);
            acb_printd(v + i, 20);
            flint_printf("\n\n");
            acb_vec_printd(v, G->phi_q, 10);
            flint_printf("\n\n");
        }
        else if (acb_rel_accuracy_bits(z) < prec - 8
                    || acb_rel_accuracy_bits(v + i) < prec - 8)
        {
                flint_printf("FAIL\n\n");
                flint_printf("q = %wu\n", q);
                flint_printf("\nL(1/2,chi_%wu(%wu,)) inaccurate\n", q, x->n);
                flint_printf("\nsingle =\n");
                acb_printd(z, 30);
                flint_printf("\ndft =\n");
                acb_printd(v + i, 30);
                flint_printf("\nerrors %ld & %ld [prec = %wu]\n",
                    acb_rel_accuracy_bits(z),
                    acb_rel_accuracy_bits(v + i), prec);
                abort();
         }

        i++;
    } while (acb_dirichlet_conrey_next(x, G) >= 0);

    acb_clear(s);
    _acb_vec_clear(v, G->phi_q);
    acb_dirichlet_char_clear(chi);
    acb_dirichlet_conrey_clear(x);
    acb_dirichlet_group_clear(G);
}

#define nq 5
#define nx 4

int main()
{

    slong i, j;
    ulong q[nq] = {3, 5, 61, 91, 800};
    ulong m[nq] = {2, 4, 11, 2, 3};
    slong prec = 150;

    acb_ptr x;
    /* cannot test at s = 1 with hurwitz */
    const char * x_r[nx] = { "0.5", "1", "1", "0.5" };
    const char * x_i[nx] = { "0", "0", "1", "6" };

    acb_t ref, res;
    /*
    default(realprecision, 100)
    X = [ 1/2, 1, 1 + I, 1/2 + 6 * I ]
    C = [Mod(2,3),Mod(4,5),Mod(11,61),Mod(2,91),Mod(3,800)]
    v = concat([ [lfun(c,x) | x<-X] | c<-C])
    apply(z->printf("\"%s\",\n",precision(real(z),70)),v)
    apply(z->printf("\"%s\",\n",precision(imag(z),70)),v)
    */
    const char * ref_r[nq * nx] = {
        "0.48086755769682862618122006323558987777682973083237186318155991033493348067023",
        "0.60459978807807261686469275254738524409468874936424685852329497846270772704212",
        "0.65552798400254803378664821634522108793943950390562746892997938695272366439543",
        "1.56831301727577320813799211138797101541772722814204756600556183229594192303169",
        "0.23175094750401575588338366176087722642788669640900589796611788733984438064542",
        "0.43040894096400403888943323295060542542457068254028965475700610399256121546113",
        "0.52127124451734699122155077366059476540647685813584432110209552194219857263213",
        "0.27554345538952180339551288674533059508689830217850843726588160302316425288119",
        "0.48926419000369574029277937487416322199001706704041739336662332096055431153634",
        "0.53487618489463947737652797478014028523937753967835874218217650866316877024685",
        "0.59822180945854055483930043368373530409360659568490328055967882287521072175030",
        "0.57333107641242898026398418236536571529256020744559201810788762987710987052691",
        "0.63562650959436738060482754500041833145501918856228134932017130492895859585857",
        "0.58238559844425623388521061809515173901502388709897277008305502466245886012927",
        "0.51027969587074040977873876733448480970861515553940454822311416220939931612454",
        "0.12930485727464247556417944278542579792607976752267116291530276949569247783180",
        "2.1717577898376043773766773888595979918343068828729776727671442492863726379365",
        "1.59480084775553338997043324643146074042666522080647790312953333780585384429160",
        "1.18088858810025653590356481638012816019876881487868657080290826197336624703523",
        "3.4156855081077462986794563990043199422106549714757808711688925562811420743252"
    };
    const char * ref_i[nq * nx] = {
        "0",
        "0",
        "0.22020604489321584265215513140893513359648656006747647412735100189927138978498",
        "-0.96945865438573217507797330416139977323695758779298609887525494651999839519954",
        "0",
        "0",
        "0.35461457326773121924283851678460530328823215076046742298881317219868652838939",
        "-0.99539202877364394787223187183283854376759830188789279576725217941822771397165",
        "-0.108902811943905225853677097712717212629591264759957601734509998247516905304233",
        "0.19365188807114177330261930603607468329782453247898970868809396136721159487113",
        "1.04370497561090171487193145841005574472705644411957863465786632716182836513300",
        "-0.23211436999860890733776901984820189055832718614668931112661490548978197438207",
        "0.0119464572932630291870372694406253796888930803905106876170081736901046822264969",
        "-0.30766179419579263291559128367018027923130717977616497873373231273162909091033",
        "-0.13330006618998077463544507824031514854466502035801914491875951675974942110581",
        "-0.56766058967929445780115371363653220980911202550251866578697696371454503377792",
        "0.97033720724583221440838051079367965353860748320561689369676918117330599741939",
        "0.38019821209845078839690749709574556766339551329040418677194351169649156337996",
        "-0.65407994257130052322391777535884554999087714891888647434188545263517159257326",
        "-1.43652482351673593824956935036654893593947145947637807266976772120632968400451"
    };

    flint_printf("l....");
    fflush(stdout);

    x = _acb_vec_init(nx);

    for (j = 0; j < nx; j++)
    {
        if (arb_set_str(acb_realref(x + j), x_r[j], prec) ||
            arb_set_str(acb_imagref(x + j), x_i[j], prec) )
        {
            flint_printf("error while setting x[%ld] <- %s+I*%s\n",
                    j, x_r[j], x_i[j]);
            abort();
        }

    }

    acb_init(ref);
    acb_init(res);

    for (i = 0; i < nq; i++)
    {
        acb_dirichlet_group_t G;
        acb_dirichlet_char_t chi;

        acb_dirichlet_group_init(G, q[i]);
        acb_dirichlet_char_init(chi, G);

        acb_dirichlet_char(chi, G, m[i]);

        for (j = 0; j < nx; j++)
        {

            if (arb_set_str(acb_realref(ref), ref_r[i * nx + j], prec) ||
                    arb_set_str(acb_imagref(ref), ref_i[i * nx + j], prec) )
            {
                flint_printf("error while setting ref <- %s+I*%s\n",
                        ref_r[i * nx + j], ref_i[i * nx + j]);
                abort();
            }
            arb_add_error_2exp_si(acb_realref(ref), -prec+10);
            arb_add_error_2exp_si(acb_imagref(ref), -prec+10);

            acb_dirichlet_l_hurwitz(res, x + j, G, chi, prec + 10);

            if (!acb_contains(ref, res))
            {
                    flint_printf("FAIL:\n\n");
                    flint_printf("q = %wu\n", q[i]);
                    flint_printf("m = %wu\n", m[i]);
                    flint_printf("x = ");
                    acb_printd(x + j, 54);
                    flint_printf("\nref = ");
                    acb_printd(ref, 54);
                    flint_printf("\nl_hurwitz(chi,x) = ");
                    acb_printd(res, 54);
                    flint_printf("\n\n");
                    abort();
            }

#if INCGAM_FIXED
            acb_dirichlet_l_incgam(res, x + j, G, chi, prec + 10);

            if (!acb_contains(ref, res))
            {
                    flint_printf("FAIL:\n\n");
                    flint_printf("q = %wu\n", q[i]);
                    flint_printf("m = %wu\n", m[i]);
                    flint_printf("x = ");
                    acb_printd(x + j, 54);
                    flint_printf("\nref = ");
                    acb_printd(ref, 54);
                    flint_printf("\nl_incgam(chi,x) = ");
                    acb_printd(res, 54);
                    flint_printf("\n\n");
                    abort();
            }
#endif

        }
        acb_dirichlet_char_clear(chi);
        acb_dirichlet_group_clear(G);

        /* test using dft */
        test_dft(q[i]);

    }
    acb_clear(ref);
    acb_clear(res);
    _acb_vec_clear(x, nx);


    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
