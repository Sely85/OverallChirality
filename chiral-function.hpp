//LQ
#include <stdio.h>
#include <cstdlib>
#include <memory>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <iterator>
#include <vector>

double* distance(double crd1[], double crd2[])
{

  double* distance = new double[3];
  for (int i=0; i<3; i++)
      distance[i]=0.0;

  distance[0] = crd2[0] - crd1[0];
  distance[1] = crd2[1] - crd1[1];
  distance[2] = crd2[2] - crd1[2];

  return distance;
}

double norm(double crd[])
{

  double norm=0.0;

  norm = sqrt(crd[0]*crd[0] + crd[1]*crd[1] + crd[2]*crd[2]);

  return norm;
}

double* crossprod(double crd1[], double crd2[])
{
  double* cross = new double[3];
  for (int i=0; i<3; i++)
      cross[i]=0.0;

  cross[0] = crd1[1]*crd2[2]-crd1[2]*crd2[1];
  cross[1] = crd1[2]*crd2[0]-crd1[0]*crd2[2];
  cross[2] = crd1[0]*crd2[1]-crd1[1]*crd2[0];

  return cross;
}

double dotprod(double crd1[], double crd2[])
{
  double dot;

  dot = crd1[0]*crd2[0] +  crd1[1]*crd2[1] + crd1[2]*crd2[2];

  return dot;
}

double single_chirality(double crd1[], double crd2[], double crd3[], double crd4[], int Na)
{
  // Calculates G0_ijkl for i,j,k,l atoms (W = weight = 1)
  //
  // G0_ijkl =  W [ ( r_ij x r_kl )r_il ]*(r_ij r_jk)(r_ij r_kl) }/ [ (|r_ij||r_jk||r_kl|)^2 |r_il| ]
  //
  // Reference: Osipov et al., Mol. Phys., 84, 1193-1206, 1995

  double*  rij = distance(crd1, crd2);
  double*  rjk = distance(crd2, crd3);
  double*  rkl = distance(crd3, crd4);
  double*  ril = distance(crd1, crd4);

  double nij = norm(rij);
  double njk = norm(rjk);
  double nkl = norm(rkl);
  double nil = norm(ril);

  double singlechi = 0.0;

  double* cross = crossprod(rij, rkl);
  double triple = dotprod(cross, ril);
  double rijjk = dotprod(rij, rjk);
  double rjkkl = dotprod(rjk, rkl);
  double num = triple*rijjk*rjkkl;

  delete[] cross;

  double den = (nij*njk*nkl)*(nij*njk*nkl)*nil;

  singlechi = num/den;

  delete[] rij;
  delete[] rjk;
  delete[] rkl;
  delete[] ril;

  return singlechi;
}

double allperm_chirality(double* crd1, double* crd2, double* crd3, double* crd4, int Na)
{
  // Compute sum of chirality index over all 24 permutations of set of 4 atoms (i, j, k, l)
  double* ri;
  double* rj;
  double* rk;
  double* rl;
  ri = new double[3];
  rj = new double[3];
  rk = new double[3];
  rl = new double[3];
  double allchi = 0.0;

  int pcnt=0;
  for (int i=0; i<4; i++)
    {
      for (int j=0; j<4; j++)
	{
	  if (i != j)
	    {
	      for (int k=0; k<4; k++)
		{
		  if (i != j && j != k && i != k)
		    {
		      for (int l=0; l<4; l++)
			{
			  if (i != j && i != k && i != l && j != k && j != l && k != l)
			    {

			      if (i == 0)
				{
				  ri[0] = crd1[0];
				  ri[1] = crd1[1];
				  ri[2] = crd1[2];
				}
			      else if (i == 1)
				{
				  rj[0] = crd1[0];
				  rj[1] = crd1[1];
				  rj[2] = crd1[2];
				}
			      else if (i == 2)
				{
				  rk[0] = crd1[0];
				  rk[1] = crd1[1];
				  rk[2] = crd1[2];
				}
			      else if (i == 3)
				{
				  rl[0] = crd1[0];
				  rl[1] = crd1[1];
				  rl[2] = crd1[2];
				}

			      if (j == 0)
				{
				  ri[0] = crd2[0];
				  ri[1] = crd2[1];
				  ri[2] = crd2[2];
				}
			      else if (j == 1)
				{
				  rj[0] = crd2[0];
				  rj[1] = crd2[1];
				  rj[2] = crd2[2];
				}
			      else if (j == 2)
				{
				  rk[0] = crd2[0];
				  rk[1] = crd2[1];
				  rk[2] = crd2[2];
				}
			      else if (j == 3)
				{
				  rl[0] = crd2[0];
				  rl[1] = crd2[1];
				  rl[2] = crd2[2];
				}

			      if (k == 0)
				{
				  ri[0] = crd3[0];
				  ri[1] = crd3[1];
				  ri[2] = crd3[2];
				}
			      else if (k == 1)
				{
				  rj[0] = crd3[0];
				  rj[1] = crd3[1];
				  rj[2] = crd3[2];
				}
			      else if (k == 2)
				{
				  rk[0] = crd3[0];
				  rk[1] = crd3[1];
				  rk[2] = crd3[2];
				}
			      else if (k == 3)
				{
				  rl[0] = crd3[0];
				  rl[1] = crd3[1];
				  rl[2] = crd3[2];
				}


			      if (l == 0)
				{
				  ri[0] = crd4[0];
				  ri[1] = crd4[1];
				  ri[2] = crd4[2];
				}
			      else if (l == 1)
				{
				  rj[0] = crd4[0];
				  rj[1] = crd4[1];
				  rj[2] = crd4[2];
				}
			      else if (l == 2)
				{
				  rk[0] = crd4[0];
				  rk[1] = crd4[1];
				  rk[2] = crd4[2];
				}
			      else if (l == 3)
				{
				  rl[0] = crd4[0];
				  rl[1] = crd4[1];
				  rl[2] = crd4[2];
				}

			      allchi = allchi + single_chirality(ri, rj, rk, rl, Na);
			      pcnt++;
			    }
			}
		    }
		}
	    }
	}
    }


  delete[] ri;
  delete[] rj;
  delete[] rk;
  delete[] rl;

  return allchi;

}
