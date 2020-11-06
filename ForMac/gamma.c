#include <stdio.h>
#include "header.h"
#include "common.h"


extern long L;
extern double Mu;

/*calculate gamma defined for approximating the conditional prob p(h | H) */
//void gammaMat(double gamma[NODE_MAX], const double THETA);
void gammaMat2(double gamma2[2][NODE_MAX], double THETA);

/*calculate the proportion of haplotype 0 at each ancestral site */
void calcProp(double *prop, long *TYPE, long ntype, long *index, long index_len);

void gammaMat(double gamma[NODE_MAX], const double THETA)
{
  short i;

  /*each element denotes the prob of mimicing mutation, when h_{x,j}=a,
  the number of ancestral types at each locus is varying*/
  for (i = 0; i < NODE_MAX; i++) {
	  gamma[i] = (double)i / ((double)i + THETA) + 0.5 * THETA / ((double)i + THETA);
  }

}

void gammaMat2(double gamma2[2][NODE_MAX], double THETA)
{
	short i;

	for (i = 0; i < NODE_MAX; i++) {
		//gamma2[0][i] = (double)i / ((double)i + THETA) + Mu * THETA / ((double)i + THETA);
		//gamma2[1][i] = (double)i / ((double)i + THETA) + (1.0-Mu) * THETA / ((double)i + THETA);
        gamma2[0][i] = (double)i / ((double)i + THETA) ;
        gamma2[1][i] = (double)i / ((double)i + THETA) ;
	}
}


void calcProp(double *prop, long *TYPE, long ntype, long *index, long index_len)
{
  /*prop is a vector of length index_len */
  long i, j, p0, p1;

  for(i=0; i<index_len; i++)
  {
    p0 = 0;
    p1 = 0;
    for(j=0; j<ntype; j++)
    {
      if(*(TYPE+(L+1)*j+index[i]) == 0) p0 += *(TYPE+(L+1)*j+L);
      if(*(TYPE+(L+1)*j+index[i]) == 1) p1 += *(TYPE+(L+1)*j+L);
    }

    if(p0+p1 > 0) prop[i] = (double)p0/((double)p0+(double)p1);
    else prop[i] = 0.5;
  }
}