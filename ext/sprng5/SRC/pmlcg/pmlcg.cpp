/*************************************************************************/
/*************************************************************************/
/*           Parallel Prime Modulus Linear Congruential Generator        */
/*                                                                       */ 
/* Modifed by: J. Ren                                                    */
/*             Florida State University                                  */
/*             Email: ren@csit.fsu.edu                                   */
/*                                                                       */
/* Based on the implementation by:                                       */
/*             Ashok Srinivasan (Apr 13, 1998)                           */
/*                                                                       */ 
/* Based on: ???                                                         */
/*                                                                       */
/*************************************************************************/
/*************************************************************************/
 
#ifdef __GNUC__
#if __GNUC__ < 4
 #include <iostream.h>
#else
 #include <iostream>
#endif /* if GNUC < 4 */
#else
  #include <iostream>
#endif /* if GNUC defined */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define NDEBUG
#include <assert.h>
#include <limits.h>
#include "memory.h"
#include "store.h"
#include "sprng.h"
#include "pmlcg.h"
#include "bignum.h"
#include "bigrat.h"
#include "basic.h"
#include <math.h>


#define MAX_STREAMS pmlcg_MAX_STREAMS
#define NGENS pmlcg_NGENS


#ifdef CONVEX
#undef _LONG_LONG  /* problems on convex compiler with 64 bit arithmetic */
#endif

#if LONG_MAX > 2147483647L 
#if LONG_MAX > 35184372088831L 
#if LONG_MAX >= 9223372036854775807L 
#define LONG_SPRNG
#define LONG64 long		/* 64 bit long */
#define store_long64 store_long
#define store_long64array store_longarray
#define load_long64 load_long
#define load_long64array load_longarray
#endif
#endif
#endif

#if !defined(LONG_SPRNG) &&  defined(_LONG_LONG)
#define LONG64 long long
#define store_long64 store_longlong
#define store_long64array store_longlongarray
#define load_long64 load_longlong
#define load_long64array load_longlongarray
#endif

#ifndef LONG64
#include "longlong.h"
#endif

#define VERSION "00"
/*** Name for Generator ***/
#define GENTYPE  VERSION "Prime modulus LCG"


#define NPARAMS 1		/*** number of valid parameters ***/
int MAX_STREAMS = (1<<30); /* Maximum number of streams for initialization */
				/* ... more streams can be spawned, though  */
   

int NGENS=0;		  /* number of random streams in current process */

#ifdef LONG64
#define limbsize 8
#else
#define limbsize 4
#endif
  
/* ************************************************************* */
/* *************************   init   ************************** */
/* ************************************************************* */

static int init(unsigned long * aa, unsigned long * x0, BigNum *k, int seed, int param)
{
  /*
     called by: initialize_int()
     calls    : 
                init_rel_prime(), prim_elt()   [ rand_lcg_mu.h ]

     params   : unsigned long aa, x0 = 'a' and 'r' arrays of a generator
                           ( empty when called )
                MP_INT k = k value to use to calculate polynomial
		param determines the power

     returns  : 'a_size' : the length of the multiplier array
	        ( also, params 'aa[]' and 'x0[]' will be filled )

        Sets up the multiplier ('a' array) and initial seed ('r' array)
        for the given value of k.   [ 'a' == 'aa' ,  'r' == 'x0' ] 
	*/

  BigNum A;
  REL_PRIME_TABLE data;
  long i, a_size;

  /* find multiplier value */
  A = (char *)MAXVAL;
  init_rel_prime(&data, &A); /* param = 2^61-1 is assumes here */
  prim_elt(&A, k, data);
  for (i=0; i<OP_SIZE; i++)
  {
    aa[i] = A.b_get_ui() & 0xffffffff;
    x0[i] = 0;
    A = b_div_2exp(A, 32);
  }

  /* initialize seed value */
  x0[0] = ((unsigned int)seed)<<1 | 1;
  free_rel_prime(&data);
  A.b_clear();

  /* calculate 'a_size' ( length of the multiplier array ) */
  i = 0;

  while (!(aa[OP_SIZE-i-1]))
    i++;

  a_size = OP_SIZE - i;

  return(a_size);
}  /* end of init() */


/* ************************************************************* */
/* *********************  initialize_int  ********************** */
/* ************************************************************* */

PMLCG **initialize(int ngen, BigNum *old_si, int seed, int param)
{
  /*
     called by: init_rng(), spawn_rng_int
     calls    : init()
     
     
     params   : int ngen = number of generators to initialize
     MP_INT old_si = value of k to use for first generator produced
     seed = encoding of starting state of generator
     param = power that determines Merssene prime
     returns  : pointer to pointers to RNGs (rngen structures)
     
     Initializes 'ngen' new generators
     ( allocates memory and gives initial values to the elements of 'rngen' )
     */

  int i,k,l,*order;
  PMLCG **q;	
  static unsigned long a[2], r[2];
  int a_size;
  
  //  order = (int *) mymalloc(ngen*sizeof(int));
  order = new int[ngen];
  /* allocate memory for 'ngen' generators */
  //  q = (PMLCG **) malloc(ngen * sizeof(PMLCG));
  q = new PMLCG * [ngen];

  if (q==NULL || order==NULL) 
    return NULL;

  for (i=0; i<ngen; i++)
  {
    q[i] = new PMLCG;

    if(q[i] == NULL)
      return NULL;

    q[i]->si = 0ul;
    q[i]->k = 0ul;
  }

  /* set up 1st generator */        
  q[0]->k = *old_si;

#ifdef LONG64
  a_size = init(a, r, &(q[0]->k),seed,param);  
  q[0]->mult = (unsigned LONG64)a[1]<<32|a[0];
  q[0]->x = (unsigned LONG64)r[1]<<32|r[0];

#else
  q[0]->a_size = init(q[0]->a, q[0]->r, &(q[0]->k),seed,param);
#endif

  q[0]->si = (*old_si) * 2;
  q[0]->si = q[0]->si + 1;

  /* set up remaining generators */
  i = 1;
  order[0] = 0;

  if (ngen>1) while (1) 
  {
    l = i;

    for (k=0; k<l; k++)
    {
      q[i]->k = q[order[k]]->si;

#ifdef LONG64
      a_size = init(a,r,&(q[i]->k),seed,param);
      q[i]->mult = (unsigned LONG64)a[1]<<32|a[0];
      q[i]->x = (unsigned LONG64)r[1]<<32|r[0];
#else
      q[i]->a_size = init(q[i]->a,q[i]->r,&(q[i]->k),seed,param);
#endif

      q[order[k]]->si = q[order[k]]->si * 2;
      q[i]->si = q[order[k]]->si;
      q[i]->si = q[i]->si + 1;

      if (ngen == ++i) 
	break;
    }

    if (ngen == i) 
      break;

    for (k=l-1; k>0; k--)
    {
      order[2*k+1] = l + k;
      order[2*k] = order[k];
    }
    order[1] = l;
  }
			
  delete [] order;
  
  return q;

} /* end of initialize */


PMLCG::PMLCG()
{
  rng_type = 5;
  gentype = NULL;
  stream_number = 0;
  nstreams = 0;
  init_seed = 0;
  parameter = 0;
  narrays = 0;
  /*** declare other variables here ***/
#ifdef LONG64
  mult = 0;
  x = 0;
#else
  r[0] = r[1] = 0;
  a[0] = a[1] = 0;
  a_size = 0;
#endif
  k = 0ul;
  si = 0ul;
}

/* Initialize random number stream */

int PMLCG::init_rng(int gn, int tg, int s, int pa)
{
/*      gives back one stream (node gennum) with updated spawning         */
/*      info; should be called total_gen times, with different value      */
/*      of gennum in [0,total_gen) each call                              */
  PMLCG **p = NULL;
  int i;
  BigNum local_k;

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
    fprintf(stderr,"ERROR - init_rng: gennum %d out of range [%d,%d).\n",
	    gn, 0, tg);
    return 0;
  }

  if (pa < 0 || pa >= NPARAMS)     /* check if parameter is valid */
  {
    fprintf(stderr,"WARNING - init_rng: parameter not valid. Using Default parameter.\n");
    pa = 0;
  }
  
  s &= 0x7fffffff;   /* Only 31 LSB of seed considered */
  local_k = Set_ui(gn); /* final seed != 0 */
  //  local_k = gn;
  p = initialize(1, &local_k, s, pa);

  if(p==NULL)
    return 0;
  else {
    rng_type = p[0]->rng_type;
    gentype = p[0]->gentype;
    stream_number = p[0]->stream_number;
    nstreams = p[0]->nstreams;
    init_seed = p[0]->init_seed;
    parameter = p[0]->parameter;
    narrays = p[0]->narrays;

#ifdef LONG64
    mult = p[0]->mult;
    x = p[0]->x;
#else
    r[0] = p[0]->r[0];
    r[1] = p[0]->r[1];
    a[0] = p[0]->a[0];
    a[1] = p[0]->a[1];
    a_size = p[0]->a_size;
#endif
    si = p[0]->si;
  }

  delete [] p;
  
  /* Initiallize data structure variables */
  rng_type = SPRNG_PMLCG;
  gentype = (char *)GENTYPE;
  stream_number = gn;
  nstreams = tg;
  init_seed = s;
  parameter = pa;
  
  narrays = 0;	/* number of arrays needed by your generator */

  while (b_cmp(si, static_cast<unsigned long int>(tg)) < 0)
    si = si * 2;

  local_k.b_clear();
  
  NGENS++;			/* NGENS = # of streams */
  
  return 1;
} 

PMLCG::~PMLCG()
{
  free_rng();
}

PMLCG::PMLCG(const PMLCG & c)
{
  rng_type = c.rng_type;
  gentype = c.gentype;
  stream_number = c.stream_number;
  nstreams = c.nstreams;
  init_seed = c.init_seed;
  parameter = c.parameter;
  narrays = c.narrays;
  /*** declare other variables here ***/
#ifdef LONG64
  mult = c.mult;
  x = c.x;
#else
  r[0] = c.r[0];
  r[1] = c.r[1];
  a[0] = c.a[0];
  a[1] = c.a[1];
  a_size = c.a_size;
#endif
  k = c.k;
  si = c.si;
}

PMLCG & PMLCG::operator= (const PMLCG & c)
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
    /*** declare other variables here ***/
#ifdef LONG64
    mult = c.mult;
    x = c.x;
#else
    r[0] = c.r[0];
    r[1] = c.r[1];
    a[0] = c.a[0];
    a[1] = c.a[1];
    a_size = c.a_size;
#endif
    k = c.k;
    si = c.si;
  }
}
/* ************************************************************* */
/* *************************  iterate  ************************* */
/* ************************************************************* */
 
void PMLCG::iterate()
{

#ifdef LONG64
  unsigned LONG64 x0, x1, x3, ul, uh, vl, vh;
#define MULT_MASK1  0x7fffffffU
#define MULT_MASK2 0x3fffffffU

#ifdef LONG_SPRNG
#define MULT_MASK3 0x1fffffffffffffffUL
#define MULT_MASK4 0x2000000000000000UL
#else
#define MULT_MASK3 0x1fffffffffffffffULL
#define MULT_MASK4 0x2000000000000000ULL
#endif

  ul = mult&MULT_MASK1;
  uh = (mult>>31)&MULT_MASK2;
  vl = x&MULT_MASK1;
  vh = (x>>31)&MULT_MASK2;
  
  x0 = ul*vl;
  x1 = ul*vh + uh*vl + (x0>>31);
  x0 &= MULT_MASK1;
  x3 = ((uh*vh)<<1) + (x1>>30);
  x0 |= (x1&MULT_MASK2)<<31;
  
  x = (x0+x3);
  
if(x&MULT_MASK4) /*Note: x != ..MASK3 since x!=0 mod prime for pmlcg */
  {
    x &= MULT_MASK3;
    x += 1;
    
    if(x == MULT_MASK4)
      x = 1;
  }
#else /* end LONG64, start 32 bit arithmetic */
/*

called by: get_rn_dbl()
calls    : add_ssaaaa(), umul_ppmm()  [ longlong.h ]

params   : int *genptr = generator to iterate
returns  : void

Performs the modular multiplication needed to iterate generator 
Xn+1 = (Xn * a) mod (2^N - 1)                                   

*/

	/*
	           aa[]    
		*  Xn[]
	       --------
               result[] -> result[] is split into 2 parts : kk[]rr[] 

	       the new Xn[] = ( kk[] + rr[] ) mod ( 2^n - 1)
	*/
  unsigned long *aa, *Xn, *rr, *kk;
  static unsigned long result[4]; /* should be atleast 2*OP_SIZE */
  static char overflow; /* should be atleast 2*OP_SIZE */
  unsigned long a0, b0, of, temp, temp2;  /* temporary storage variables */
  unsigned long prod_lo, prod_hi, res_lo, res_hi;
  long i,j;  /* counter variables           */
  int param = 0;
  
  aa = a;
  Xn = r;

  memset(result,0,4*sizeof(unsigned long)); /* initialize to 0 */
  overflow = 0;
  

  /* result[] = aa[] * Xn[] */
  a0 = aa[0];
  b0 = Xn[0];
  umul_ppmm(prod_hi,prod_lo, a0,b0); 
  result[0] = prod_lo;
  result[1] = prod_hi;
    
  b0 = Xn[1];
  umul_ppmm(prod_hi,prod_lo, a0,b0); 
  res_lo = result[1];
  add_ssaaaa(of,temp2, 0,prod_lo, 0,res_lo);
  result[1] = temp2;
  add_ssaaaa(res_hi,res_lo,0,prod_hi, 0,of);
  result[2] = res_lo;
  overflow = res_hi;

  if (a_size == 2)
  {
    a0 = aa[1];
    b0 = Xn[0];
    res_lo = result[1];
    res_hi = result[2];
    umul_ppmm(prod_hi,prod_lo, a0,b0); 
    add_ssaaaa(of,temp2, 0,prod_lo, 0,res_lo);
    result[1] = temp2;
    add_ssaaaa(temp,temp2,  0, prod_hi,  0,res_hi);
    add_ssaaaa(res_hi,res_lo,  temp,temp2, 0,of);
    result[2] = res_lo;
    overflow += res_hi;

    b0 = Xn[1];
    res_lo = result[2];
    res_hi = result[3];
    umul_ppmm(prod_hi,prod_lo, a0,b0); 
    add_ssaaaa(of,temp2, 0,prod_lo, 0,res_lo);
    result[2] = temp2;
    add_ssaaaa(temp,of, 0,of, 0,overflow);
    add_ssaaaa(temp,temp2,  0, prod_hi,  0,res_hi);
    add_ssaaaa(res_hi,res_lo,  temp,temp2, 0,of);
    result[3] = res_lo;
  }
    
  /*  rr = low(result) (R)    kk = hi(result) (K)  */
  rr = result;
  kk = result + OP_SIZE;

  /* shift 'kk' left */
  temp2 = 0;
  for (i=0; i<OP_SIZE; i++)
  {
    temp = kk[i];
    kk[i] = ((temp<<RNGBITS)&0xffffffff) + temp2;
    temp2 = temp >> SHIFT;
  }

  /* move extra bits at top of rr[] into start of kk[] */
  temp = rr[OP_SIZE-1];
  temp2 = temp >> SHIFT;
  rr[OP_SIZE-1] = ((temp<<RNGBITS)&0xffffffff) >> RNGBITS;
  kk[0] += temp2;

  /*  Xn+1 = rr + kk   */
  a0 = 0;

  for (i=0; i<OP_SIZE; i++)
  {
    temp = rr[i];
    temp2 = kk[i];
    add_ssaaaa(of,b0, 0,temp, 0,temp2);
    add_ssaaaa(temp2,temp, of,b0,  0,a0);
    Xn[i] = temp;
    a0 = temp2;
  }

  /*  perform mod operation  Xn+1 = Xn+1 mod 2^n - 1 */
  /*  Xn+1 = ( r & (2^n-1) ) + ( r >> n )  */
  temp2 = Xn[OP_SIZE - 1] >> SHIFT;
  Xn[OP_SIZE - 1] &= MASK; 

  for (i=0; i<OP_SIZE; i++)
  {
    temp = Xn[i]; 
    add_ssaaaa(temp2,res_lo, 0,temp, 0,temp2);
    Xn[i] = res_lo;
  }
#endif
  
} /* end of iterate() */


double PMLCG::get_rn_dbl()
{
#ifdef LONG64
  static double dtemp[1] = {0.0};

#ifdef LONG_SPRNG
#define EXPO 0x3ff0000000000000UL
#else
#define EXPO 0x3ff0000000000000ULL
#endif

  iterate();

#if defined(CONVEX) || defined(O2K) || defined(SGI) || defined(GENERIC)
  *((unsigned LONG64 *) dtemp) = (genptr->x>>9) | EXPO;
  return *dtemp - (double) 1.0;
#else
  return (x>>9)*2.2204460492503131e-16;
#endif

#else  /* 32 bit arithmetic */
  double num1,num2;
  long i;

  iterate();

  num1 = (double) r[0];
  num2 = (double) r[1];
  num2 *= (double) 0XFFFFFFFF + 1.0;
  num1 += num2;	
	
  num2 = (double) 0XFFFFFFFF + 1.0;
  num2 *= (double) 0X1FFFFFFF + 1.0;
  num2 -= 1.0;
  num1 /= num2;        
  
  return (num1);
#endif
  
} 



/* Return a random integer */

int PMLCG::get_rn_int()
{
#ifdef LONG64
  iterate();
  return (int) (x>>30);
#else
  unsigned long irn;
  
  iterate();

  irn = (r[1]<<2) | ((r[0]&0xc0000000)>>30);
  
  return (int) (irn&0x7fffffff); 	
#endif
  
} 


/* Return a single precision random number */

float PMLCG::get_rn_flt()
{
  return (float) get_rn_dbl();
}


/*************************************************************************/
/*************************************************************************/
/*                  SPAWN_RNG: spawns new generators                     */
/*************************************************************************/
/*************************************************************************/

int PMLCG::spawn_rng(int nspawned, Sprng ***newgens)
{
  PMLCG ** genptr;
  int i;
  
  if (nspawned <= 0) /* is nspawned valid ? */
  {
    nspawned = 1;
    fprintf(stderr,"WARNING - spawn_rng: nspawned <= 0. Default value of 1 used for nspawned\n");
  }
  
  genptr = initialize(nspawned, &si,init_seed,parameter);

  if(genptr == NULL)	   /* allocate memory for pointers to structures */
  {
    *newgens = NULL;
    return 0;
  }
  else
  {
    *newgens = (Sprng **) genptr;

    for(i=0; i<nspawned; i++)
    {
      genptr[i]->rng_type = rng_type;
      genptr[i]->gentype = (char *)GENTYPE;
      genptr[i]->stream_number = stream_number;
      genptr[i]->nstreams = nstreams;
      genptr[i]->init_seed = init_seed; 
      genptr[i]->parameter = parameter;
  
      genptr[i]->narrays = 0;		/* number of arrays needed by your generator */

      NGENS++;
    }
  }
  
  return nspawned;
}


/* Free memory allocated for data structure associated with stream */

int PMLCG::free_rng()
{
  int i;
  
  assert(this != NULL);
  
  k.b_clear();
  si.b_clear();

  //  free(this);

  NGENS--;
  return NGENS;
}


int PMLCG::pack_rng(char **buffer)
{

  int i, size;
  unsigned char *p, *initp;
  unsigned int m[2], temp;
  
  /*std::cout << "packing the following:\n";
  std::cout << init_seed << " = seed " << parameter << " = parameter " << stream_number << " = streamnum\n";
  std::cout << mult << " = mult " << x << " = x\n\n";
  std::cout << "k = " << k << "\n";
  std::cout << "si = " << si << "\n\n";*/

  #ifdef LONG64
  size = 4 + 3*4 + 2*8;// strlen(gentype)+1;
  #else
  size = 4 + 3*4 + 2*8 + 4;
  #endif
  
  size += 2 * sizeof(si.size) + si.size * sizeof(long int) + k.size * sizeof(long int) + 2;
  /* The new load/store routines make using sizeof unnecessary. Infact, */
  /* using sizeof could be erroneous. */
  initp = p = (unsigned char *) mymalloc(size);

  if(p == NULL)
  {
    *buffer = NULL;
    return 0;
  }
  
  p += store_int(rng_type,4,p);

  p += store_int(init_seed,4,p);
  p += store_int(stream_number,4,p);
  p += store_int(parameter,4,p);
  
  #ifdef LONG64
  p += store_long64(mult,8,p);
  p += store_long64(x,8,p);
  #else
  p += store_int(a_size,4,p);
  m[0] = r[0]>>8;
  m[1] = r[0]<<24 | r[1];
  p += store_intarray(m,2,4,p);
  
  m[0] = a[0]>>8;
  m[1] = a[0]<<24 | a[1];
  p += store_intarray(m,2,4,p);

  #endif
  
  /* Storing bignum member vars */
  p += store_long(si.size, sizeof(si.size), p);
  p += store_longarray(si.v, si.size, sizeof(long int), p);
  p += 1;
  memcpy (p, &si.sign, 1);
  
  
  p += store_long(k.size, sizeof(k.size), p);
  p += store_longarray(k.v, k.size, sizeof(long int), p);
  p += 1;
  memcpy (p, &k.sign, 1);
  
  *buffer = (char *) initp;
  assert(p-initp == size);
  
  return p-initp;
  
  /*char *temp_buffer;
  int size, i;
  int pos=0;
  size = 4 + sizeof(PMLCG) + narrays * sizeof(int) + strlen(gentype) + 1;
  size += k.size * limbsize;
  size += si.size * limbsize;
  
  temp_buffer = new char[size];

  if(temp_buffer == NULL)
  {
    *buffer = NULL;
    return 0;
  }
  
  pos += store_int(rng_type,4,(unsigned char*)(temp_buffer+pos));
  /*strcpy(temp_buffer+pos,gentype);
  pos += strlen(gentype)+1;
  */
  
  
  /*memcpy(temp_buffer+pos,this,sizeof(PMLCG));
  pos += sizeof(PMLCG);
  
  memcpy(temp_buffer+pos, k.v, k.size * limbsize);
  pos += k.size * limbsize;
  memcpy(temp_buffer+pos, si.v, si.size * limbsize);
  pos += si.size * limbsize;

  assert(pos == size);
  
  *buffer = temp_buffer;
  */
  
  
  return size;
}


int PMLCG::unpack_rng(char *packed)
{
  unsigned int i, m[2];
  unsigned char *p;
  
  p = (unsigned char *) packed;

  if(this == NULL) 
    return 0;
	
  p += load_int(p,4,(unsigned int *)&rng_type);
  
  if (rng_type != SPRNG_PMLCG) {
    fprintf(stderr,"ERROR: Unpacked ' %d ' instead of ' %d '\n", rng_type, SPRNG_LCG64); 
    return 0;
  }
  p += load_int(p,4,(unsigned int *)&init_seed);
  p += load_int(p,4,(unsigned int *)&stream_number);
  p += load_int(p,4,(unsigned int *)&parameter);


#ifdef LONG64                   /* 64 bit integer types */
  p += load_long64(p,8,&mult);
  p += load_long64(p,8,&x);
#else  /* No 64 bit type available */
  p += load_int(p,4,&a_size);

  p += load_intarray(p,2,4,&m);
  r[1] = m[1]&0xffffff; 
  r[0] = m[1]>>24 | m[0]<<8;
  
  p += load_intarray(p,2,4,&m);  
  a[1] = m[1]&0xffffff; 
  a[0] = m[1]>>24 | m[0]<<8;
#endif

  /* Loading bignum member vars */
  p += load_long(p, sizeof(si.size), &(si.size));
  p += load_longarray(p, si.size, sizeof(long int), si.v);
  p += 1;
  memcpy (&si.sign, p, 1);

  p += load_long(p, sizeof(k.size), &(k.size));
  p += load_longarray(p, k.size, sizeof(long int), k.v);
  p += 1;
  memcpy (&k.sign, p, 1);


  narrays = 0;  
    
  NGENS++;

  /*std::cout << "unpacked the following:\n";
  std::cout << std::dec << init_seed << " = seed " << parameter << " = parameter " << stream_number << " = streamnum\n";
  std::cout << mult << " = mult " << x << " = x\n\n";
  std::cout << "k = " << k << "\n";
  std::cout << "si = " << si << "\n\n";*/
  return 1;


/*
  int i;
  int pos=0;

  if (this == NULL)
    return 0;

 // pos += 4; // skip rng_type 
  /*if(strcmp(packed+pos,GENTYPE) != 0)
  {
    fprintf(stderr,"ERROR: Unpacked ' %.24s ' instead of ' %s '\n",  
	    packed+pos, GENTYPE); 
    return 0;
  }
  else
    gentype = (char *)GENTYPE;

  pos += strlen(gentype)+1;*/
  /*
  pos += load_int(p,4,(unsigned int *) &rng_type)
  if (rng_type != SPRNG_PMLCG) {
    fprintf(stderr,"ERROR: Unpacked ' %d ' instead of ' %d '\n", id, SPRNG_LCG); 
    return 0;
  }
    
  memcpy(this,packed+pos,sizeof(PMLCG));
  pos += sizeof(PMLCG);
    
  //  q->k._mp_d = (mp_limb_t *) mymalloc(q->k._mp_alloc*sizeof(mp_limb_t));
  //  q->si._mp_d = (mp_limb_t *) mymalloc(q->si._mp_alloc*sizeof(mp_limb_t));

  k.v = new unsigned long int[k.size];
  si.v = new unsigned long int[si.size];

  //  if(q->k._mp_d == NULL || q->si._mp_d == NULL)
  if (k.v == NULL || si.v == NULL)
    return 0;

  /*
  memcpy(q->k._mp_d,packed+pos,q->k._mp_alloc*sizeof(mp_limb_t));
  pos += q->k._mp_alloc*sizeof(mp_limb_t);
  memcpy(q->si._mp_d,packed+pos,q->si._mp_alloc*sizeof(mp_limb_t));
  pos += q->si._mp_alloc*sizeof(mp_limb_t);
  */
/*
  memcpy(k.v, packed+pos, k.Size() * limbsize);
  pos += k.Size() * limbsize;
  memcpy(si.v, packed+pos, si.Size() * limbsize);
  pos += si.Size() * limbsize;

  NGENS++;

  return 1;
  */
}

      
int PMLCG::get_seed_rng()
{
  return init_seed;
}


int PMLCG::print_rng()
{
  printf("\n%s\n", GENTYPE+2);
  printf("\n \tseed = %d, stream_number = %d\tparameter = %d\n\n", init_seed, stream_number, parameter);	
  /*#ifdef LONG64
  printf("multiplier = %llu, seed = %llu\n", mult, x);
  #endif*/

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
