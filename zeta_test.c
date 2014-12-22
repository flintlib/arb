
#include "/home/therickaman/Documents/SagePrograms/arb/arb.h"
#include "/home/therickaman/Documents/SagePrograms/arb/acb.h"
#include "/home/therickaman/Documents/SagePrograms/arb/acb_poly.h"
#include "/home/therickaman/Documents/SagePrograms/arb/fmpr.h"
#include "mpfr.h"

struct hurwitz_params;

typedef int (*acb_func)(const struct hurwitz_params *context,acb_poly_t output,acb_poly_t input);

struct hurwitz_params{  //Should be easily generalizable
  int init;             //Is everything initialized(that needs to be)?
  acb_t a;              //
  long int precision;   //
  long int derivatives; //
  acb_func funct;       //Function to be integrated
};                      //

void hurwitz_params_clear(struct hurwitz_params *context){
	acb_clear(context->a);
	flint_free(context);
}

int lambda_hurwitz(const struct hurwitz_params *context,acb_poly_t out,acb_poly_t inp_s){
  acb_poly_zeta_series(out,inp_s,context->a,0,context->derivatives,context->precision);
  return 0;
}

struct hurwitz_params *build_hurwitz_struct(acb_t inp_a,long int inp_d,long int inp_prec){

  struct hurwitz_params *result = flint_malloc(sizeof(struct hurwitz_params));
  
  if(result->init == 0){
        acb_init(result->a);
  }


  acb_set(result->a,inp_a);
  result->precision=inp_prec;
  result->derivatives=inp_d;
  result->funct = lambda_hurwitz;
  result->init = 1;
  return  result;
}

int main(int argc, char *argv[]) {

  int i;
  long int deriv = 10;
  long int prec  = 1000;

  acb_t a;
  acb_t z;
  acb_poly_t s;
  acb_poly_t out;

  acb_init(a);
  acb_init(z);
  acb_poly_init(s);
  acb_poly_init(out);

  acb_set_ui(a,(long int) 1);
  acb_div_ui(a,a,(long int) 2,prec);
  struct hurwitz_params *p = build_hurwitz_struct(a,deriv,prec);
  acb_set_ui(z,2);

  acb_poly_set_coeff_si(s, 1, 1);
  p->funct(p,out,s);
  acb_poly_printd(out,20);
  printf("\n");

  acb_poly_set_acb(s,z);
  acb_poly_set_coeff_si(s, 1, 1);
  p->funct(p,out,s);
  acb_poly_printd(out,20);
  printf("\n");

  acb_clear(a);
  acb_clear(z);
  acb_poly_clear(s);
  acb_poly_clear(out);
  return 1;
}
