/*************************************************************************/
/*************************************************************************/
/*           Parallel Multiplicative Lagged Fibonacci Generator          */
/*                                                                       */ 
/* Modifed by: J. Ren                                                    */
/*             Florida State University                                  */
/*             Email: ren@csit.fsu.edu                                   */
/*                                                                       */
/* Based on the implementation by:                                       */
/*             Ashok Srinivasan (1997)                                   */
/*************************************************************************/
/*************************************************************************/

#if __GNUC__ > 3
 #include <iostream>
#else
 #include <iostream.h>
#endif
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define NDEBUG
#include <assert.h>
#include "memory.h"
#include "sprng.h"
#include "mlfg.h"
#include "store.h"

#define MAX_STREAMS mlfg_MAX_STREAMS
#define NGENS mlfg_NGENS
#define valid mlfg_valid

#define VERSION "00"
/*** Name for Generator ***/
#define GENTYPE VERSION "Multiplicative Lagged Fibonacci Generator"

#define NPARAMS 11		/* number of valid parameters */
int MAX_STREAMS = 0x7fffffff;  /* Maximum number streams to init_sprng */

struct vstruct {
      int L;
      int K;
      int LSBS;     /* number of least significant bits that are 1 */
      int first;    /* the first seed whose LSB is 1 */
};

const struct vstruct valid[] = { {17,5,1,10}, {31,6,1,2},
{55,24,1,11}, {63,31,1,14}, {127,97,1,21}, {521,353,1,100},
{521,168,1,83}, {607,334,1,166}, {607,273,1,105}, {1279,418,1,208}, {1279,861,1,233} };


/* int64.h (START) */

#if LONG_MAX > 2147483647L 
#if LONG_MAX > 35184372088831L 
#if LONG_MAX >= 9223372036854775807L 
#define LONG_SPRNG
#define LONG64 long             /* 64 bit long */
#endif
#endif
#endif

#if !defined(LONG_SPRNG) &&  defined(_LONG_LONG)
#define LONG64 long long
#endif


#ifdef LONG64

typedef unsigned LONG64 uint64;

#ifdef LONG_SPRNG
#define MINUS1 0xffffffffffffffffUL             /* -1 (mod 2^(BITS-2)) */
#define ONE    0x1UL
#define MASK64 0xffffffffffffffffUL
#else
#define MINUS1 0xffffffffffffffffULL            /* -1 (mod 2^(BITS-2)) */
#define ONE    0x1ULL
#define MASK64 0xffffffffffffffffULL
#endif

#define multiply(a,b,c) {c = (a)*(b); c &= MASK64;}
#define add(a,b,c)      {c = (a)+(b); c &= MASK64;}
#define decrement(a,c)  {c = (a)-1; c &= MASK64;}
#define And(a,b,c)      {c = (a)&(b);}
#define Or(a,b,c)       {c = (a)|(b);}
#define Xor(a,b,c)      {c = (a)^(b);}
#define notzero(a)      (a==0?0:1)
#define lshift(a,b,c)   {c = (a)<<(b); c &= MASK64;} /* b is an int */
#define rshift(a,b,c)   {c = (a)>>(b); c &= MASK64;} /* b is an int */
#define highword(a)     ((unsigned int)((a)>>32))
#define lowword(a)      ((unsigned int)((a)&0xffffffff))
#define set(a,b)        {b = (a)&MASK64;}
#define seti(a,b)       {b = (a)&MASK64;} /* b is an int */
#define seti2(a,b,c)    {c = (b); c <<= 32; c |= (a); c &= MASK64;}/*a,b=+int*/

#else   /* Simulate 64 bit arithmetic on 32 bit integers */

typedef unsigned int uint64[2];

static const uint64 MASK64={0xffffffffU,0xffffffffU};
static const uint64 MINUS1={0xffffffffU,0xffffffffU}; /* -1 (mod 2^(BITS-2)) */
static uint64 ONE={0x1U,0x0U}; 
#define TWO_M32 2.3283064365386962e-10 /* 2^(-32) */

#define And(a,b,c)      {c[0] = a[0]&b[0]; c[1] = a[1]&b[1];}
#define Or(a,b,c)       {c[0] = a[0]|b[0]; c[1] = a[1]|b[1];}
#define Xor(a,b,c)      {c[0] = a[0]^b[0]; c[1] = a[1]^b[1];}
#define notzero(a)      ((a[0]==0 && a[1]==0)?0:1)
#define multiply(a,b,c) {c[1] = a[0]*b[1]+a[1]*b[0];\
                           c[1] += (unsigned int) (((double)a[0]*(double)b[0])\
                                                                *TWO_M32);\
                           c[0] = a[0]*b[0]; And(c,MASK64,c);}
#define add(a,b,c)      {unsigned int t = a[0]+b[0]; \
                           c[1] = a[1]+b[1]+(t<a[0]);c[0]=t;\
                           And(c,MASK64,c);}
#define decrement(a,c)  {if(a[0]==0){c[1]=a[1]-1;c[0]=0xffffffff;} \
                            else c[0] = a[0]-1; And(c,MASK64,c);}

static void lshift(uint64 a,int b,uint64 c)   
{
  if(b<32)
    {c[1] = (a[1]<<b)|(a[0]>>(32-b)); c[0] = a[0]<<(b);}
  else {c[1]=a[0]<<(b-32);c[0]=0;} 
  And(c,MASK64,c);
} 

static void rshift(uint64 a,int b,uint64 c)   
{
  if(b<32)
    {c[0] = (a[0]>>b)|(a[1]<<(32-b));c[1] = a[1]>>(b);}
  else {c[0]=a[1]>>(b-32);c[1]=0;} 
  And(c,MASK64,c);
} 

#define highword(a)     ((a)[1])
#define lowword(a)      ((a)[0])
#define set(a,b)        {b[0] = a[0];b[1]=a[1];And(b,MASK64,b);}
#define seti(a,b)       {b[0] = a;b[1]=0;} /* b is an int */
#define seti2(a,b,c)    {c[1] = b; c[0] = a; And(c,MASK64,c);}/*a,b = +ve int*/

#endif  /* LONG64 or 32 bit */  


static int store_uint64(uint64 l, unsigned char *c)
{
  int i;
  unsigned int m[2];
  
  m[0] = highword(l);
  m[1] = lowword(l);
  
  c += store_intarray(m,2,4,c);

  return 8;             /* return number of chars filled */
}


static int store_uint64array(uint64 *l, int n, unsigned char *c)
{
  int i;
  
  for(i=0; i<n; i++)
    c += store_uint64(l[i],c);

  return 8*n;
}


int load_uint64(unsigned char *c, uint64 *l)
{
  int i;
  unsigned int m[2];
  
  c += load_intarray(c,2,4,m);
  seti2(m[1],m[0],(*l));
  
  return 8;
}


int load_uint64array(unsigned char *c, int n, uint64 *l)
{
  int i;
  
  for(i=0; i<n; i++)
    c += load_uint64(c,l+i);

  return 8*n;
}

/* int64.h (END) */


#define TWO_M52 2.2204460492503131e-16 /* 2^(-52) */
#define TWO_M64 5.4210108624275222e-20 /* 2^(-64) */
#define BITS 62			/* Initialization of ALFG part is m-2 bits */
#define MAX_BIT_INT (BITS-2)
#define RUNUP (2*BITS)		/* Do RUNUP iterations after initialization */
#define GS0 0x372f05ac

#ifdef LONG64
#define INT_MOD_MASK (MASK64>>(64-BITS))
#define INT_MASK (MASK64>>(64-BITS+1))
#define INTX2_MASK ((((uint64)1)<<MAX_BIT_INT)-1)
static const uint64 SEED_MASK=0x5a38; 
//static uint64 SEED_MASK=0x5a38; 
#else
static const uint64 INT_MOD_MASK={0xffffffffU,0x3fffffffU};
static const uint64 INT_MASK={0xffffffffU,0x1fffffffU};
static const uint64 INTX2_MASK={0xffffffffU,0x0fffffffU};
static const uint64 SEED_MASK={0x5a38,0x0U};
/*
static uint64 INT_MOD_MASK={0xffffffffU,0x3fffffffU};
static uint64 INT_MASK={0xffffffffU,0x1fffffffU};
static uint64 INTX2_MASK={0xffffffffU,0x0fffffffU};
static uint64 SEED_MASK={0x5a38,0x0U};
*/
#endif


int NGENS=0;		  /* number of random streams in current process */


/*************************************************************************/
/*************************************************************************/
/*            ROUTINES USED TO CREATE GENERATOR FILLS                    */
/*************************************************************************/
/*************************************************************************/

static int bitcnt( int x)
{
  unsigned i=0,y;

  for (y=(unsigned)x; y; y &= (y-1) ) 
    i++;

  return(i);
}

static void advance_reg(int *reg_fill)
{
/*      the register steps according to the primitive polynomial         */
/*           (64,4,3,1,0); each call steps register 64 times             */
/*      we use two words to represent the register to allow for integer  */
/*           size of 32 bits                                             */

  const int mask = 0x1b;
  int adv_64[4][2];
  int i,new_fill[2];
  unsigned temp;

  adv_64[0][0] = 0xb0000000;
  adv_64[0][1] = 0x1b;
  adv_64[1][0] = 0x60000000;
  adv_64[1][1] = 0x2d;
  adv_64[2][0] = 0xc0000000;
  adv_64[2][1] = 0x5a;
  adv_64[3][0] = 0x80000000;
  adv_64[3][1] = 0xaf;
  new_fill[1] = new_fill[0] = 0;
  temp = mask<<27;

  for (i=27;i>=0;i--) 
  {
    new_fill[0] = (new_fill[0]<<1) | (1&bitcnt(reg_fill[0]&temp));
    new_fill[1] = (new_fill[1]<<1) | (1&bitcnt(reg_fill[1]&temp));
    temp >>= 1;
  }

  for (i=28;i<32;i++) 
  {
    temp = bitcnt(reg_fill[0]&(mask<<i));
    temp ^= bitcnt(reg_fill[1]&(mask>>(32-i)));
    new_fill[0] |= (1&temp)<<i;
    temp = bitcnt(reg_fill[0]&adv_64[i-28][0]);
    temp ^= bitcnt(reg_fill[1]&adv_64[i-28][1]);
    new_fill[1] |= (1&temp)<<i;
  }

  reg_fill[0] = new_fill[0];
  reg_fill[1] = new_fill[1];
}

static void get_fill(uint64 *n, uint64 *r, int param, unsigned seed)
{
  /*int i,j,k,temp[2], length;*/
  unsigned int i,j,k,temp[2], length;
  uint64 tempui;
  
  length = valid[param].L;
  
  /* Initialize the shift register with the node number XORed with seed    */
  temp[1] = highword(n[0]);
  temp[0] = lowword(n[0])^seed;
  if (!temp[0])
    temp[0] = GS0;


  //  advance_reg(temp); /* Advance the shift register some */
  //  advance_reg(temp);
  int temp_int[2];
  temp_int[0] = static_cast<int>(temp[0]);
  temp_int[1] = static_cast<int>(temp[1]);
  advance_reg(temp_int);
  advance_reg(temp_int);

  /* The first word of the RNG is defined by the LSBs of the node number   */
  And(n[0],INT_MASK,tempui);
  lshift(tempui,1,r[0]);
  
  /* The RNG is filled with the bits of the shift register, at each time   */
  /* shifted up to make room for the bits defining the canonical form;     */
  /* the node number is XORed into the fill to make the generators unique  */

  for (i=1;i<length-1;i++) 
  {
    //    advance_reg(temp);
    advance_reg(temp_int);
  
    seti2(temp[0],temp[1],tempui);
    Xor(tempui,n[i],tempui);
    And(tempui,INT_MASK,tempui);
    lshift(tempui,1,r[i]);
  }
  seti(0,r[length-1]);
/*      the canonical form for the LSB is instituted here                */
  k = valid[param].first + valid[param].LSBS;

  for (j=valid[param].first;j<k;j++)
    Or(r[j],ONE,r[j]);

  return;
}


/* left shift array 'b' by one, and place result in array 'a' */

static void si_double(uint64 *a,  uint64 *b, int length)
{
  int i;
  uint64 mask1, temp1;
  
  lshift(ONE,MAX_BIT_INT,mask1);
  
  And(b[length-2],mask1,temp1);

  if (notzero(temp1))
    fprintf(stderr,"WARNING: si_double -- RNG has branched maximum number of times.\n\t Independence of generators no longer guaranteed\n");

  And(b[length-2],INTX2_MASK,temp1);
  lshift(temp1,1,a[length-2]);

  for (i=length-3;i>=0;i--) 
  {
    And(b[i],mask1,temp1);

    if(notzero(temp1)) 
      add(a[i+1],ONE,a[i+1]);

    And(b[i],INTX2_MASK,temp1);
    lshift(temp1,1,a[i]);
  }
}

static void pow3(uint64 n, uint64 *ui)		/* return 3^n (mod 2^BITS) */
{
  uint64 p, value, temp, bit, temp2, temp3;
  int exponent;
  
  set(n,p);
  seti(3,temp);
  seti(1,temp3);
  And(n,temp3,temp2);
  
  if (notzero(temp2)) {
    seti(3,value);
  }
  else {
    seti(1,value);
  }

  seti(1,bit);
  
  for(exponent=2; exponent<64; exponent++)
  {
    multiply(temp,temp,temp);
    lshift(bit,1,bit);
    
    And(bit,n,temp2);
    if(notzero(temp2))
      multiply(value,temp,value);
  }
  
  And(value,MASK64,value);
  
  set(value,(*ui));
}

static void findseed(int sign, uint64 n, uint64 *ui)
{
  uint64 temp;
  
  pow3(n,&temp);
  
  if(sign&1)
    multiply(temp,MINUS1,temp);
  
  set(temp,(*ui));
}

void MLFG::advance_state()
{
  int lv = lval, kv = kval;
  int lptr;
  hptr--;
  
  if(hptr < 0)
    hptr = lv-1;

  lptr = hptr + kv;

  if (lptr>=lv)
   lptr -= lv;

  multiply(lags[hptr],lags[lptr],lags[hptr]);
}

static MLFG** initialize(int rng_type_local, int ngen_local, int param_local, 
			 unsigned int seed_local, uint64 *nstart_local, 
			 unsigned int initseed_local)
{
  int i,j,k,l,m,*order, length;
  MLFG ** q;
  uint64 *nindex, temp1, mask;

  length = valid[param_local].L;
  
  //  order = (int *) mymalloc(ngen_local*sizeof(int));
  //  q = (MLFG **) mymalloc(ngen_local*sizeof(MLFG *));
  order = new int[ngen_local];
  q = new MLFG * [ngen_local];

  if (q == NULL || order == NULL) 
    return NULL;

  for (i=0;i<ngen_local;i++) 
  {
    q[i] = new MLFG;

    if (q[i] == NULL) 
      return NULL;

    q[i]->rng_type = rng_type_local;
    q[i]->hptr = 0;   /* This is reset to lval-1 before first iteration */
    //    q[i]->si = (uint64 *) mymalloc((length-1)*sizeof(uint64));
    //    q[i]->lags = (uint64 *) mymalloc(length*sizeof(uint64));
    q[i]->si = new uint64[length-1];
    q[i]->lags = new uint64[length];
    q[i]->lval = length;
    q[i]->kval = valid[param_local].K;
    q[i]->parameter = param_local;
    q[i]->seed = seed_local;
    q[i]->init_seed = initseed_local;
    q[i]->narrays=2;
    q[i]->gentype = (char *)GENTYPE;
    
    if (q[i]->lags == NULL || q[i]->si == NULL) 
      return NULL;
  }
/*      specify register fills and node number arrays                    */
/*      do fills in tree fashion so that all fills branch from index     */
/*           contained in nstart array                                   */
  q[0]->stream_number = lowword(nstart_local[0]);
  get_fill(nstart_local,q[0]->lags,param_local,seed_local);
  si_double(q[0]->si,nstart_local,length);

  set(ONE,mask);

  for(m=0; m<length; m++)
  {
    And(SEED_MASK,mask,temp1);

    if(notzero(temp1))
      findseed(1,q[0]->lags[m], &q[0]->lags[m]);
    else
      findseed(0,q[0]->lags[m], &q[0]->lags[m]);

    lshift(mask,1,mask);
  }
  
  add(q[0]->si[0],ONE,q[0]->si[0]);

  i = 1;
  order[0] = 0;

  if (ngen_local>1) 
    while (1) 
    {
      l = i;
      for (k=0;k<l;k++) 
      {
	nindex = q[order[k]]->si;
	q[i]->stream_number = lowword(nindex[0]);
	get_fill(nindex,q[i]->lags,param_local,seed_local);
	si_double(nindex,nindex, length);
	for (j=0;j<length-1;j++) 
	  set(nindex[j],q[i]->si[j]);

	set(ONE,mask);

	for(m=0; m<length; m++)
	{
	  And(SEED_MASK,mask,temp1);

	  if(notzero(temp1))
	    findseed(1,q[i]->lags[m], &q[i]->lags[m]);
	  else
	    findseed(0,q[i]->lags[m], &q[i]->lags[m]);

	  lshift(mask,1,mask);
	}
	
	add(q[i]->si[0],ONE,q[i]->si[0]);

	if (ngen_local == ++i) 
	  break;
      }
      
      if (ngen_local == i) 
	break;
                
      for (k=l-1;k>0;k--) 
      {
	order[2*k+1] = l+k;
	order[2*k] = order[k];
      }
      order[1] = l;
    }

  free(order);

  for (i=ngen_local-1;i>=0;i--) 
  {
    k = 0;
    for (j=1;j<length-1;j++)
      if (notzero(q[i]->si[j])) 
	k = 1;
    if (!k) 
      break;
    for (j=0;j<length*RUNUP;j++)
      q[i]->advance_state();
  }

  while (i>=0)
  {
    for (j=0;j<4*length;j++)
      q[i]->advance_state();
    i--;
  }   

  return q;
}

/* Initialize random number stream */

MLFG::MLFG()
{
  rng_type = SPRNG_MLFG;
  gentype = NULL;
  stream_number = 0;
  nstreams = 0;
  init_seed = 0;
  parameter = 0;
  narrays = 0;
  array_sizes = NULL;
  arrays = NULL;
  lags = NULL;
  si = NULL;
  hptr = 0;
  lval = 0;
  kval = 0;
  seed = 0;
}

int MLFG::init_rng(int gn, int tg, int s, int pa)
{
/*      gives back one stream (node gennum) with updated spawning         */
/*      info; should be called total_gen times, with different value      */
/*      of gennum in [0,total_gen) each call                              */
  MLFG **p=NULL;
  uint64 *nstart=NULL,*si_local;
  int i, length, k;
  
  if (tg <= 0) /* Is total_gen valid ? */
  {
    tg = 1;
    fprintf(stderr,"WARNING - init_rng: Total_gen <= 0. Default value of 1 used for total_gen\n");
  }

  if (gn >= MAX_STREAMS) /* check if gen_num is valid    */
    fprintf(stderr,"WARNING - init_rng: gennum: %d > maximum number of independent streams: %d\n\tIndependence of streams cannot be guranteed.\n",
	    gn, MAX_STREAMS); 

  if (gn < 0 || gn >= tg) /* check if gen_num is valid    */
  {
    fprintf(stderr,"ERROR - init_rng: gennum %d out of range [%d,%d).\n", gn, 0, tg); 
    return 0;
  }

  s &= 0x7fffffff;		/* Only 31 LSB of seed considered */

  if (pa < 0 || pa >= NPARAMS)     /* check if parameter is valid */
  {
    fprintf(stderr,"WARNING - init_rng: parameter not valid. Using Default parameter.\n");
    pa = 0;
  }
  
  length = valid[pa].L; /* determine parameters   */
  k = valid[pa].K;
  
  /*      define the starting vector for the initial node                  */
  //`   nstart = (uint64 *) mymalloc((length-1)*sizeof(uint64));
  nstart = new uint64[length-1];

  if (nstart == NULL)
    return 0;

  seti(gn,nstart[0]);

  for (i=1;i<length-1;i++) 
    seti(0,nstart[i]);

  p = initialize(4, 1,pa,s^GS0,nstart,s);  /* create a generator  */

  if (p==NULL) 
    return 0;
  else
  {
    rng_type = p[0]->rng_type;
    gentype = p[0]->gentype;
    stream_number = p[0]->stream_number;
    nstreams = p[0]->nstreams;
    init_seed = p[0]->init_seed;
    parameter = p[0]->parameter;
    narrays = p[0]->narrays;
    arrays = p[0]->arrays;
    lags = p[0]->lags;
    si = p[0]->si;
    hptr = p[0]->hptr;
    lval = p[0]->lval;
    kval = p[0]->kval;
    seed = p[0]->seed;

    free(p);
  }
  
/*      update si array to allow for future spawning of generators       */
  si_local = si;

  while (lowword(si_local[0]) < tg && !highword(si_local[0])) 
    si_double(si_local,si_local,length);

  NGENS++;
  
  rng_type = 4;
  stream_number = gn;
  nstreams = tg;
  
  //  array_sizes = (int *) mymalloc(narrays*sizeof(int));
  //  arrays = (int **) mymalloc(narrays*sizeof(int *));
  array_sizes = new int[narrays];
  arrays = new int * [narrays];

  if(array_sizes == NULL || arrays == NULL)
    return 0;

  arrays[0] = (int *) lags;
  arrays[1] = (int *) si;
  array_sizes[0] = lval*sizeof(uint64)/sizeof(int);
  array_sizes[1] = (lval-1)*sizeof(uint64)/sizeof(int);

  return 1;
} 

MLFG::~MLFG()
{
  free_rng();
}

MLFG::MLFG(const MLFG & c)
{
  rng_type = c.rng_type;
  gentype = c.gentype;
  stream_number = c.stream_number;
  nstreams = c.nstreams;
  init_seed = c.init_seed;
  parameter = c.parameter;
  narrays = c.narrays;
  array_sizes = c.array_sizes;
  arrays = c.arrays;
  lags = c.lags;
  si = c.si;
  hptr = c.hptr;
  lval = c.lval;
  kval = c.kval;
  seed = c.seed;
}

MLFG & MLFG::operator= (const MLFG & c)
{
  if (this != &c) {
    this->free_rng();

    rng_type = c.rng_type;
    gentype = c.gentype;
    stream_number = c.stream_number;
    nstreams = c.nstreams;
    init_seed = c.init_seed;
    parameter = c.parameter;
    narrays = c.narrays;
    array_sizes = c.array_sizes;
    arrays = c.arrays;
    lags = c.lags;
    si = c.si;
    hptr = c.hptr;
    lval = c.lval;
    kval = c.kval;
    seed = c.seed;
  }
}

/* Returns a double precision random number */

double MLFG::get_rn_dbl()
{
  advance_state();	
  /*printf("\t sprng: %lu\n", lags[hptr]);*/
#ifdef LONG64  
  return (lags[hptr]>>12)*TWO_M52;
#else
  return (lags[hptr][1])*TWO_M32 + (lags[hptr][0])*TWO_M64;
#endif
} 

/* Return a random integer */

int MLFG::get_rn_int()
{
  advance_state();

#ifdef LONG64
  return lags[hptr]>>33;
#else
  return lags[hptr][1]>>1;
#endif
} 

/* Return a single precision random number */

float MLFG::get_rn_flt()
{
  /* If you have a more efficient way of computing the random integer,
     then please replace the statement below with your scheme.        */

  return (float) get_rn_dbl();
}


/*************************************************************************/
/*************************************************************************/
/*                  SPAWN_RNG: spawns new generators                     */
/*************************************************************************/
/*************************************************************************/

int MLFG::spawn_rng(int nspawned, Sprng ***newgens)
{
  MLFG ** genptr;
  int i;
  uint64 *p;
  
  if (nspawned <= 0) /* is nspawned valid ? */
  {
    nspawned = 1;
    fprintf(stderr,"WARNING - spawn_rng: nspawned <= 0. Default value of 1 used for nspawned\n");
  }
  
  p = si;

  genptr = initialize(rng_type, nspawned,parameter,seed,p,init_seed);

  if(genptr == NULL)	   /* allocate memory for pointers to structures */
  {
    *newgens = NULL;
    return 0;
  }
  
  si_double(p,p,lval);

  for(i=0; i<nspawned; i++)
  {
    //    genptr[i]->array_sizes = (int *) mymalloc(genptr[i]->narrays*sizeof(int));
    //    genptr[i]->arrays = (int **) mymalloc(genptr[i]->narrays*sizeof(int *));
    genptr[i]->array_sizes = new int[genptr[i]->narrays];
    genptr[i]->arrays = new int * [genptr[i]->narrays];

    if(genptr[i]->array_sizes == NULL || genptr[i]->arrays == NULL)
      return 0;

    genptr[i]->arrays[0] = (int *)(genptr[i]->lags);
    genptr[i]->arrays[1] = (int *)(genptr[i]->si);
    genptr[i]->array_sizes[0] = genptr[i]->lval*sizeof(uint64)/sizeof(int);
    genptr[i]->array_sizes[1] = (genptr[i]->lval-1)*sizeof(uint64)/sizeof(int);
  }
  
  NGENS += nspawned;
      
  *newgens = (Sprng **) genptr;

  return nspawned;
}


/* Free memory allocated for data structure associated with stream */

int MLFG::free_rng()
{
  int i;
  
  assert(this != NULL);

  for(i=0; i<narrays; i++)
    delete [] arrays[i];

  if(narrays > 0)
  {
    delete [] array_sizes;
    delete [] arrays;
  }

  //  free(this);  
  NGENS--;

  return NGENS;
}



int MLFG::pack_rng(char **buffer)
{
  unsigned char *p, *initp;
  int size, i;

  size = 4 + 24+16*lval + strlen(gentype)+1;
  
  initp = p = (unsigned char *) mymalloc(size); /* allocate memory */

  if(p == NULL)
  {
    *buffer = NULL;
    return 0;
  }
  
  p += store_int(rng_type,4,p);
  /*strcpy((char *) p,gentype);
  p += strlen(gentype)+1;*/
  //p += store_int(SPRNG_MLFG,4,p);
  p += store_int(stream_number,4,p);
  p += store_int(nstreams,4,p);
  p += store_int(init_seed,4,p);
  p += store_int(parameter,4,p);
  p += store_int(hptr,4,p);
  p += store_int(lval,4,p);
  p += store_int(kval,4,p);
  p += store_int(seed,4,p);
  p += store_uint64array(si,lval-1,p);
  p += store_uint64array(lags,lval,p);
  
  *buffer = (char *) initp;
  assert(p-initp == size);
  
  return p-initp;
}


int MLFG::unpack_rng(char *packed)
{
  unsigned char *p;
  int i, found;

  if(this == NULL) 
    return 0;

  p = (unsigned char *) packed;
  p += load_int(p,4,(unsigned int *)&rng_type);

  /*if(strcmp((char *) p,GENTYPE) != 0)
  {
    fprintf(stderr,"ERROR: Unpacked ' %.24s ' instead of ' %s '\n", p, GENTYPE); 
    return 0;
  }
  else
    gentype = (char *)GENTYPE;

  p += strlen(gentype)+1;    */

  if (rng_type != SPRNG_MLFG) {
    fprintf(stderr,"ERROR: Unpacked ' %d ' instead of ' %d '\n", rng_type, SPRNG_MLFG); 
    return 0;
  }
  p += load_int(p,4,(unsigned int *)&stream_number);
  p += load_int(p,4,(unsigned int *)&nstreams);
  p += load_int(p,4,(unsigned int *)&init_seed);
  p += load_int(p,4,(unsigned int *)&parameter);
  p += load_int(p,4,(unsigned int *)&hptr);
  p += load_int(p,4,(unsigned int *)&lval);
  p += load_int(p,4,(unsigned int *)&kval);
  p += load_int(p,4,(unsigned int *)&seed);

/*      check values of parameters for consistency                       */
  for(i=found=0; i<NPARAMS; i++)
    if(lval==valid[i].L && kval==valid[i].K)
    {
      found = 1;
      break;
    }
  
  if(found == 0)
  {
    fprintf(stderr,"ERROR: Unpacked parameters are not acceptable.\n");
    exit(EXIT_FAILURE);
  }

  narrays = 2;
  //  array_sizes = (int *) mymalloc(narrays*sizeof(int));
  //  arrays = (int **) mymalloc(narrays*sizeof(int *));
  array_sizes = new int[narrays];
  arrays = new int * [narrays];

  if(array_sizes == NULL || arrays == NULL) 
    return 0;

  array_sizes[0] = lval*sizeof(uint64)/sizeof(int);
  array_sizes[1] = (lval-1)*sizeof(uint64)/sizeof(int);
  
  for(i=0; i<narrays; i++)
  {
    //    arrays[i] = (int *) mymalloc(array_sizes[i]*sizeof(int));
    arrays[i] = new int[array_sizes[i]];

    if(arrays[i] == NULL) 
      return 0;
  }   

  lags = (uint64 *) arrays[0];
  si   = (uint64 *) arrays[1];
  p += load_uint64array(p,lval-1,si);
  p += load_uint64array(p,lval,lags);  
    
  NGENS++;

  return 1;
}
      

int MLFG::get_seed_rng()
{
  return init_seed;
}

int MLFG::print_rng()
{
  printf("\n%s\n", GENTYPE+2);
  printf("\n \tseed = %d, stream_number = %d\tparameter = %d\n\n", init_seed, stream_number, parameter);

  return 1;
}




/***********************************************************************************
* SPRNG (c) 2014 by Florida State University                                       *
*                                                                                  *
* SPRNG is licensed under a                                                        *
* Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. *
*                                                                                  *
* You should have received a copy of the license along with this                   *
* work. If not, see <http://creativecommons.org/licenses/by-nc-sa/4.0/>.           *
************************************************************************************/
