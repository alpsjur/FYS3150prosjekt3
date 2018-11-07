//#include "mpi.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <armadillo>

#include "metropolis.hpp"

ofstream ofile;
void writeToFile(double *values, double T, int mcs, int count);
void solveForDifferentMCS(double T, bool ordered, string filename);
void calculateProbability(double T, string filename);

int main(int argc, char const *argv[]) {
  /*
  string file1 = "../data/EMvsMCS_ordered_1.dat";
  string file2 = "../data/EMvsMCS_unordered_1.dat";
  string file3 = "../data/EMvsMCS_ordered_24.dat";
  string file4 = "../data/EMvsMCS_unordered_24.dat";


  solveForDifferentMCS(1, true, file1);
  solveForDifferentMCS(1, false, file2);
  solveForDifferentMCS(2.4, true, file3);
  solveForDifferentMCS(2.4, false, file4);
  */

  string file5 = "../data/probability1.dat";
  string file6 = "../data/probability24.dat";

  calculateProbability(1, file5);
  //calculateProbability(2.4, file6);

  return 0;
}

void writeToFile(double *values, double T, int mcs, int count){
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) << mcs;
  ofile << setw(15) << setprecision(8) << T;
  ofile << setw(15) << setprecision(8) << values[0];
  ofile << setw(15) << setprecision(8) << values[1];
  ofile << setw(15) << setprecision(8) << values[2];
  ofile << setw(15) << setprecision(8) << values[3];
  ofile << setw(15) << setprecision(8) << values[4];
  ofile << setw(15) << setprecision(8) << count << endl;
}

void solveForDifferentMCS(double T, bool ordered, string filename){
  int L = 20; double p[1000];
  ofile.open(filename);

  double J = 1; double k = 1; double beta = 1/(k*T);
  double values[5]; long idum = -1; int count;

  for (int mcs = 1e0; mcs <= 1e6; mcs*=10){
    solveGivenT(L, mcs, T, k, J, values, idum, ordered, count,p);
    calculateVarNormalize(values, L, mcs);
    calculateCChi(values, k, T);
    writeToFile(values, T, mcs, count);
  }
  ofile.close();
}

void calculateProbability(double T, string filename){
  int L = 20; int mcs = 1e6; bool ordered = false;
  double J = 1; double k = 1; double beta = 1/(k*T);
  double values[5]; long idum = -1; int count;
  double p[1000] = {0};

  ofile.open(filename);

  solveGivenT(L, mcs, T, k, J, values, idum, ordered, count,p);
  calculateVarNormalize(values, L, mcs);
  ofile << values[0] << ' ' << values[1] << endl;
  for (int i = 0; i<1000;++i){
    if (p[i] != 0){
      ofile << setw(5) << setprecision(3) << -i/(double)(L*L);
      ofile << setw(15) << setprecision(8) << p[i] << endl;
    }
  }
  ofile.close();
}
