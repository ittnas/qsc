#define MAXVAL "152D02C7E14AF6800000" /*maximum # of independent streams*/

struct BIGNUM_ARRAY_TYPE;

#define  OP_SIZE    2
#define  RNGBITS       3
#define  SHIFT      29
#define  MASK       0X1FFFFFFF
#define  MAGIC_NUM  "222222222222222"
#define  MAGIC_DEN  "60454F554C0000"
#define  PRIM       37
#define  POWER_N    61

/* Values pertain to this particular parameter: 2^61-1 as modulus*/
struct BIGNUM_ARRAY_TYPE init_factors();


/***********************************************************************************
* SPRNG (c) 2014 by Florida State University                                       *
*                                                                                  *
* SPRNG is licensed under a                                                        *
* Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. *
*                                                                                  *
* You should have received a copy of the license along with this                   *
* work. If not, see <http://creativecommons.org/licenses/by-nc-sa/4.0/>.           *
************************************************************************************/
