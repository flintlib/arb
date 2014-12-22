#include "/home/therickaman/Documents/SagePrograms/arb/arb.h"
#include "acb.h"
#include "acb_poly.h"
#include "fmpr.h"
#include "mpfr.h"

// Compile With:
// gcc -o zeta_values zeta_values.c -lgmp -lmpfr -lflint -lflint -larb -lm -I/usr/include/flint

void acb_printVals(const acb_t z, long digits){

	arf_printd(arb_midref(acb_realref(z)), digits);
	if (arf_sgn(arb_midref(acb_imagref(z))) < 0)
	{
		arf_t t;
		arf_init_neg_shallow(t, arb_midref(acb_imagref(z)));
		printf(" - ");
		arf_printd(t, digits);
	}
	else
	{
		printf(" + ");
		arf_printd(arb_midref(acb_imagref(z)), digits);
	}
	printf("*I");
}

int getInput(int *argc,char *argv[],acb_t *s,acb_t *a,long int *numDerivs,long int *prec){

    	if (*argc < 10)
    	{
		printf("Please Call This Function With:\n");
        	printf("zeta_values [-real x] [-imag y] [-a a] [-numDerivs numDerivs] [-prec prec]\n");
        	return 1;
    	}

	//These cases are easy.  prec and numDerivs are always integer.
	//Also, we need prec to know the precision of x,y,a.

	if(!strcmp(argv[9],"-prec"))
	    sscanf(argv[10],"%lu",prec);

	if(!strcmp(argv[7],"-numDerivs"))
	    sscanf(argv[8],"%lu",numDerivs);


	mpfr_t x,y,mpfr_a;
	arf_t real,imag,arf_a;

	mpfr_init2(x,(int) *prec);
	mpfr_init2(y,(int) *prec);
	mpfr_init2(mpfr_a,(int) *prec);
	arf_init(real);
	arf_init(imag);
	arf_init(arf_a);

	if (!strcmp(argv[1], "-real"))
	    mpfr_set_str (x,argv[2],10,MPFR_RNDN);

	if(!strcmp(argv[3],"-imag"))
	    mpfr_set_str (y,argv[4],10,MPFR_RNDN);
	
	if(!strcmp(argv[5],"-a"))
	    mpfr_set_str (mpfr_a,argv[6],10,MPFR_RNDN);

	arf_set_mpfr(real,x);
	arf_set_mpfr(imag,y);
	arf_set_mpfr(arf_a,mpfr_a);

	arb_set_arf(acb_realref(*s),real);
	arb_set_arf(acb_imagref(*s),imag);
	arb_set_arf(acb_realref(*a),arf_a);

	mpfr_clear(x);
	mpfr_clear(y);
	mpfr_clear(mpfr_a);
	arf_clear(real);
	arf_clear(imag);
	arf_clear(arf_a);
}

void gaussQuad(acb_ptr *weights[],acb_ptr *nodes[],int points){

	//TODO: DO THIS FOR GENERAL QUADRATURE RULE!

	arb_set_d(weights[0],(double) 0.4179591836734694);
	arb_set_d(weights[1],(double) 0.3818300505051189);
	arb_set_d(weights[2],(double) 0.3818300505051189);
	arb_set_d(weights[3],(double) 0.2797053914892766);
	arb_set_d(weights[4],(double) 0.2797053914892766);
	arb_set_d(weights[5],(double) 0.1294849661688697);
	arb_set_d(weights[6],(double) 0.1294849661688697);

	arb_set_d(nodes[0],(double) 0.0);
	arb_set_d(nodes[1],(double) 0.4058451513773972);
	arb_set_d(nodes[2],(double) -0.4058451513773972);
	arb_set_d(nodes[3],(double) -0.7415311855993945);
	arb_set_d(nodes[4],(double) 0.7415311855993945);
	arb_set_d(nodes[5],(double) -0.9491079123427585);
	arb_set_d(nodes[6],(double) 0.9491079123427585);
}

(*arb_calc_func_t)(arb_ptr out,
    const arb_t inp, void * param, long order, long prec);

int sinx(acb_ptr out, const acb_t inp, void * params, long order, long prec)
{
    int xlen = FLINT_MIN(2, order);
    acb_set(out, inp);
    if (xlen > 1)
        acb_one(out + 1);
    _acb_poly_sin_series(out, out, xlen, order, prec);
    return 0;
}

int zetaDeriv(arb_ptr out,const arb_t inp, void * param, long order, long prec){

}


void findNumZerosInBox(arb_t x1,arb_t x2,arb_t T1,arb_t T2,void (*func)(acb_t)){

	//This algorithm takes as input *acb_t and returns the value though the parameter.


	acb_t x;
	acb_init(x);
	acb_set_ui(x,1);

	void f(*acb_t) = (*func)(*acb_t);

	f(&x);

	acb_printd(x,10);





}

int zeroRegions(acb_t a,int deriv,int minT,int maxT){

	

	int j,loopCond;
	long int numDerivs=(long int)deriv+1;//THIS ENSURES THAT THE 0,1,...,deriv,deriv+1 derivatives are computed.
	long int prec=300,digits = prec*0.301029995663981-1,sizeSeven = 7;

	arb_ptr x;
	arb_ptr w;

	arb_ptr param;

	acb_ptr s_left;
	acb_ptr s_bottom;
	acb_ptr s_right;
	acb_ptr s_top;
	acb_ptr left;
	acb_ptr bottom;
	acb_ptr right;
	acb_ptr top;
	
	arb_t x1;
	arb_t x2;
	arb_t T1;
	arb_t T2;
	arb_t maxX;
	arb_t minX;
	arb_t h1;
	arb_t c0;
	
	arb_init(x1);
	arb_init(x2);
	arb_init(T1);
	arb_init(T2);
	arb_init(maxX);
	arb_init(minX);
	arb_init(h1);
	arb_init(c0);

	x        = _arb_vec_init(sizeSeven);
	w        = _arb_vec_init(sizeSeven);
	param    = _arb_vec_init(sizeSeven);

	s_left   = _acb_vec_init(sizeSeven);
	s_bottom = _acb_vec_init(sizeSeven);
	s_right  = _acb_vec_init(sizeSeven);
	s_top    = _acb_vec_init(sizeSeven);
	left     = _acb_vec_init(sizeSeven);
	bottom   = _acb_vec_init(sizeSeven);
	right    = _acb_vec_init(sizeSeven);
	top      = _acb_vec_init(sizeSeven);

	arb_set_d(maxX,(double) 1.13588*numDerivs+2.0);
	arb_set_d(h1,(double) 0.33534);
	arb_neg(minX,maxX);

	arb_set_d(T1,(double) -.1);
	arb_set_d(T2,(double) .1);
	arb_set(x1,minX);
	arb_add(x2,x1,h1,prec);

	arb_sub(c0,T1,T2,prec);   //c0 = T1-T2
	arb_div_ui(c0,c0,2,prec); //c0 = (T1-T2)/2


	findNumZerosInBox(T1,T2,T1,T2,void (*testfunc));

	return 1;

	loopCond = (arf_cmp(arb_midref(x1),arb_midref(x2))<0)*(arb_is_negative(x2)!=0);

	acb_vec_set_im_vector(s_left,x,c0,sizeSeven);//Set sleft=((T1-T2)/2)*x[i]*I   //
	//acb_vec_set(s_right,s_left,sizeSeven); //Set sright = sleft = ((T1-T2)/2)*x[i]*I

	while(loopCond){
		arb_mul_si(x2,x2,(long int) arb_is_negative(x2)!=0,prec);

		acb_vec_set_re_arb(s_left,x1,sizeSeven); //Set s_left = x1 + ((T1-T2)/2)*x[i]*I
		acb_vec_set_re_arb(s_right,x2,sizeSeven);//Set s_right = x2 + ((T1-T2)/2)*x[i]*I


		//arb_set_arf(acb_realref(s),real);
		//arb_set_arf(acb_imagref(s),imag);

		//left


		//printf("%lf,%lf\n",x1,x2);
		//arb_printd(x1,20);
		//printf(",");
		//arb_printd(x2,20);
		//printf("\n");
		
		

		arb_set(x1,x2);
		arb_add(x2,x2,h1,prec);
		loopCond = arb_is_negative(x1)!=0;
		//TODO: CHANGE C0 when x2 = 0.

	}	
	



	/////////////////////////////////////
	//CLEAR ALL USED VARIABLES FROM ARB//
	/////////////////////////////////////

	//acb_clear(s);
	//_acb_vec_clear(zetaDeriv,numDerivs);
	//_acb_vec_clear(zetaRatios1,numDerivs-1);
	//_acb_vec_clear(zetaRatios2,numDerivs-1);

}

int main(int argc, char *argv[])
{

	long int i,digits,numDerivs,prec;

	acb_t s;
	acb_t a;
   	acb_ptr result;

	acb_init(s);
	acb_init(a);
	
	getInput(&argc,argv,&s,&a,&numDerivs,&prec);

	numDerivs++;

	digits = (long int) prec*0.301029995663981;

	result = _acb_vec_init(numDerivs);

	_acb_poly_zeta_cpx_series(result,s,a,0,numDerivs,prec);

	zeroRegions(a,100,-1000,1000);


	//for (i = 0; i < 7; i++)
    	//{
	     //arb_printd(w+i,16);
	     //printVals(result + i,digits); 
	     //printf("\n");
	     //arb_printd(n+i,16);
	     //printf("\n");
    	//}

	acb_clear(s);
	acb_clear(a);
	_acb_vec_clear(result, numDerivs);

    	return 0;
}
