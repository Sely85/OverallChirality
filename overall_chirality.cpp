#include <stdio.h>
#include <cstdlib>
#include <memory>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include "chiral-function.hpp"
 
using namespace std;

//Define a class with specie, x, y, z
class coordinates
{
public:
  coordinates(std::string sp, double xxx, double yyy, double zzz): specie(sp),x(xxx), y(yyy), z(zzz)
  {}
  coordinates(): specie(), x(0.0), y(0.0), z(0.0)
  {}
  std::string specie;
  double x;
  double y;
  double z;
};

int main(int argc, char *argv[])
{
  unsigned const stacksize = 500;

  /***************************************  
             Read input file
  ***************************************/
  std::cout << " " << "\n";
  std::cout << "Welcome to the OverallChirality code!"  << "\n";
  std::cout << "   Usage: " << "\n";
  std::cout << "          ./OverallChirality <xyz> " << "\n";
  std::cout << " " << "\n";
  std::cout << " " << "\n";

  ifstream infile (argv[1]);

  if ( !infile.is_open() )
    {
      std::cout << "ERROR: Could not open file" << "\n";
      return 0;
    }

  if ( argc < 2) 
    {
      std::cout << "ERROR: Missing input parameter" << "\n";
      return 0;
    }


  string null;
  int totatom;

  infile >> totatom;
  infile >> null;

  std::cout << " Particles......" << totatom << "\n";

  std::vector<coordinates> coord(totatom);
  //read xyz input data
  for (int i=0; i<totatom; i++)
    {
      infile >> coord[i].specie >> coord[i].x >> coord[i].y >> coord[i].z;
    }    

  int Na = totatom;
  double overall = 0.0;

  int cnt = 0;

#pragma omp parallel
  for (int i=0; i<Na-3; i++)
    {
      for (int j=i+1; j<Na-2; j++)
	{
	  for (int k=j+1; k<Na-1; k++)
	    {

              #pragma omp for reduction(+ : overall) 
	      for (int l=k+1; l<Na; l++)
		{
		  double* r1 = new double[3];
		  double* r2 = new double[3];
		  double* r3 = new double[3];
		  double* r4 = new double[3];
		  r1[0] = coord[i].x;
		  r1[1] = coord[i].y;
		  r1[2] = coord[i].z;
		  r2[0] = coord[j].x;
		  r2[1] = coord[j].y;
		  r2[2] = coord[j].z;
		  r3[0] = coord[k].x;
		  r3[1] = coord[k].y;
		  r3[2] = coord[k].z;
		  r4[0] = coord[l].x;
		  r4[1] = coord[l].y;
		  r4[2] = coord[l].z;
		  overall = overall + allperm_chirality(r1, r2, r3, r4, Na);
		  cnt++;
		  delete[] r1;
		  delete[] r2;
		  delete[] r3;
		  delete[] r4;
		}
	    }
	}
    }
  std::cout << " " << "\n";
  std::cout << " Overall chirality....... "  << overall/totatom  << "\n";
  std::cout << " " << "\n";
} 


