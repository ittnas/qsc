#include <cstdio>
#include <cstdlib>
#include "cputime.h"

using namespace std;

#define TIMING_TRIAL_SIZE 1000000

int main()
{
  int i;
  double temp1, temp2, rn;
  double temp_mult = TIMING_TRIAL_SIZE/1.0e6;
  
  temp1 = cputime();

  for(i=0; i<TIMING_TRIAL_SIZE; i++)
    rn = drand48();
  
  temp2 = cputime();
  

  if(temp2-temp1 < 1.0e-15 )
  {
    printf("Timing Information not available/not accurate enough.\n\t...Exiting\n");
    exit(1);
  }
  
  /* The next line is just used to ensure that the optimization does not remove the call to the RNG. Nothing is really printed.             */
  fprintf(stderr,"Last random number generated\n", rn);
  
  printf("\nUser + System time Information (Note: MRS = Million Random Numbers Per Second)\n");
  printf("\tDRAND48:\t Time = %7.3f seconds => %8.4f MRS\n", 
	 temp2-temp1, temp_mult/(temp2-temp1));
  putchar('\n');
  
  return 0;
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
