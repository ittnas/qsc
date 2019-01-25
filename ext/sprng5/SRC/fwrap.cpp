/************************************************************************/
/************************************************************************/
/*                                                                      */
/*      This package of C++ wrappers is intended to be called from a    */
/*      FORTRAN program. The main purpose of the package is to mediate  */
/*      between the call-by-address and call-by-value conventions in    */
/*      the two languages. In most cases, the arguments of the C++      */
/*      routines and the wrappers are the same.                         */
/*                                                                      */
/*      The wrappers should be treated as FORTRAN function calls.       */
/*      Their names will be the same, except in FORTRAN they will       */
/*      not have the final _F77. (Example: ffree_rng_F77 becomes        */
/*      ffree_rng                                                       */
/*                                                                      */
/* Note: This code is a header file to facilitte inlining.              */
/************************************************************************/
/************************************************************************/

extern "C"
{
  #include "fwrap.h"
}

#include "sprng_cpp.h"
#include "memory.h"

#if __GNUC__ > 3
 #include <iostream>
#else
 #include <iostream.h>
#endif
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>


extern "C"
{
  int fget_seed_rng_F77(int **genptr);
  int ffree_rng_F77(int **genptr);
  int fmake_new_seed_F77(void);
  int * finit_rng_sim_F77(int *seed,  int *mult, int *gtype = 0);
  int * finit_rng_F77(int * rng_type, int *gennum, int *total_gen, int *seed, int *length);
  int fspawn_rng_F77(int **genptr,  int *nspawned, int **newGen);
  int fget_rn_int_sim_F77(void);
  int fget_rn_int_F77(int **genptr);
  float fget_rn_flt_sim_F77(void);
  float fget_rn_flt_F77(int **genptr);
  double fget_rn_dbl_sim_F77(void);
  double fget_rn_dbl_F77(int **genptr);
  int fpack_rng_F77(int **genptr, char *buffer);
  int fpack_rng_simple_F77(char *buffer);
  int * funpack_rng_F77(char *buffer, int *rng_type);
  int * funpack_rng_simple_F77(char *buffer, int *rng_type);
  int fprint_rng_F77( int **genptr);
  int fprint_rng_simple_F77(void);
  int fseed_mpi_F77(void);

  int * finit_rng_simmpi_F77(int *seed, int *mult, int * rng_type);
  int fget_rn_int_simmpi_F77(void);
  float fget_rn_flt_simmpi_F77(void);
  double fget_rn_dbl_simmpi_F77(void);

  double fcpu_t_F77(void);
}


double fcpu_t_F77(void)
{
  double   current_time;

#ifdef RUSAGE_SELF
  struct rusage temp;

  getrusage(RUSAGE_SELF, &temp);

  current_time = (temp.ru_utime.tv_sec + temp.ru_stime.tv_sec +
		  1.0e-6*(temp.ru_utime.tv_usec + temp.ru_stime.tv_usec));

#elif defined(CLOCKS_PER_SEC)
  current_time = clock()/((double) CLOCKS_PER_SEC);

#else
  fprintf(stderr,"\nERROR: Timing routines not available\n\n");
  current_time = 0.0;
#endif

  return (current_time);
}

int fget_seed_rng_F77(int **genptr)
{
  Sprng * ptr = (Sprng *) *genptr;
  
  return ptr->get_seed_rng();
}

int ffree_rng_F77(int **genptr)
{
  Sprng * ptr = (Sprng *) *genptr;
  
  return ptr->free_rng();
}

int fmake_new_seed_F77(void)
{
  return make_new_seed();
}

int * finit_rng_sim_F77(int *seed,  int *mult, int *gtype)
{
  return init_rng_simple(*seed, *mult, *gtype);
}

int * finit_rng_F77(int * rng_type, int *gennum, int *total_gen, int *seed, int *length)
{
  Sprng * ptr = SelectType(* rng_type);
  int status = ptr->init_rng(*gennum, *total_gen, *seed, *length);
  
  if (status == 1)
    return (int *) ptr;
  else
    return (int *) NULL;
}

int fspawn_rng_F77(int **genptr,  int *nspawned, int **newGen)
{
  int i, n;
  Sprng ** tmpGen;
  Sprng * ptr = (Sprng *) *genptr;
  
  n =  ptr->spawn_rng(*nspawned, &tmpGen);
  
  for (i=0; i< n; i++) 
    newGen[i] = (int *) tmpGen[i];

  if(n != 0)
    delete [] tmpGen;
  
  return n;
}

int fget_rn_int_sim_F77(void)
{
  return get_rn_int_simple();
}

int fget_rn_int_F77(int **genptr)
{
  Sprng * ptr = (Sprng *) *genptr;

  return ptr->get_rn_int();
}

float fget_rn_flt_sim_F77(void)
{
  return get_rn_flt_simple();
}

float fget_rn_flt_F77(int **genptr)
{
  Sprng * ptr = (Sprng *) *genptr;
  
  return ptr->get_rn_flt();
}

double fget_rn_dbl_sim_F77(void)
{
  return get_rn_dbl_simple();
}

double fget_rn_dbl_F77(int **genptr)
{
  Sprng * ptr = (Sprng *) *genptr;

  return ptr->get_rn_dbl();
}

int fpack_rng_F77(int **genptr, char *buffer)
{
  int size;
  char *temp;
  Sprng * ptr = (Sprng *) *genptr;

  size = ptr->pack_rng(&temp);

  if(temp != NULL)
    {
      memcpy(buffer,temp,size);
      delete [] temp;
    }
  
  return size;
}

int fpack_rng_simple_F77(char *buffer)
{
  int size;
  char *temp;

  size = pack_rng_simple(&temp);

  if(temp != NULL)
    {
      memcpy(buffer,temp,size);
      delete [] temp;
    }
  
  return size;
}

int * funpack_rng_F77(char *buffer, int * rng_type)
{
  Sprng * ptr = SelectType(* rng_type);
  int status = ptr->unpack_rng(buffer);

  if (status == 1)
    return (int *) ptr;
  else
    return (int *) NULL;
}

int * funpack_rng_simple_F77(char *buffer, int * rng_type)
{
  return unpack_rng_simple(buffer, *rng_type);
}

int fprint_rng_F77( int **genptr)
{
  Sprng * ptr = (Sprng *) *genptr;

  return ptr->print_rng();
}

int fprint_rng_simple_F77(void)
{
  return print_rng_simple();
}




int fseed_mpi_F77(void)
{
#ifdef SPRNG_MPI
  return make_new_seed_mpi();
#else
  return -1;
#endif
}



int * finit_rng_simmpi_F77(int *seed, int *mult, int *rng_type)
{
  return init_rng_simple_mpi(*seed, *mult, *rng_type);
}

int fget_rn_int_simmpi_F77(void)
{
  return get_rn_int_simple_mpi();
}

float fget_rn_flt_simmpi_F77(void)
{
  return get_rn_flt_simple_mpi();
}

double fget_rn_dbl_simmpi_F77(void)
{
  return get_rn_dbl_simple_mpi();
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
