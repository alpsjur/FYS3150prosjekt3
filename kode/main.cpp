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
void solveForDifferentMCS(double T, bool ordered, string filename, int, int, int);
void solveForDifferentMCS2(double T, bool ordered, string filename, int, int, int);
void calculateProbability(double T, string filename);

int main(int argc, char const *argv[]) {


  string file1 = "../data/EMvsMCS_ordered_1.bin";
  string file2 = "../data/EMvsMCS_unordered_1.bin";
  string file3 = "../data/EMvsMCS_ordered_24.bin";
  string file4 = "../data/EMvsMCS_unordered_24.bin";


  int mcsStart=100, mcsStop=3e4, mcsStep=100;
  solveForDifferentMCS(1, true, file1, mcsStart, mcsStop, mcsStep);
  solveForDifferentMCS(1, false, file2, mcsStart, mcsStop, mcsStep);
  solveForDifferentMCS(2.4, true, file3, mcsStart, mcsStop, mcsStep);
  solveForDifferentMCS(2.4, false, file4, mcsStart, mcsStop, mcsStep);

/*

  string file5 = "../data/count_ordered_1.bin";
  string file6 = "../data/count_unordered_1.bin";
  string file7 = "../data/count_ordered_24.bin";
  string file8 = "../data/count_unordered_24.bin";

  int mcsStart=10, mcsStop=1e6, mcsStep=10;
  solveForDifferentMCS2(1, true, file5, mcsStart, mcsStop, mcsStep);
  solveForDifferentMCS2(1, false, file6, mcsStart, mcsStop, mcsStep);
  solveForDifferentMCS2(2.4, true, file7, mcsStart, mcsStop, mcsStep);
  solveForDifferentMCS2(2.4, false, file8, mcsStart, mcsStop, mcsStep);


  string file9 = "../data/probability1.dat";
  string file10 = "../data/probability24.dat";

  calculateProbability(1, file9);
  calculateProbability(2.4, file10);
*/
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

void solveForDifferentMCS(double T, bool ordered, string filename, int mcsStart, int mcsStop, int mcsStep){
  int L = 20; double p[1000]; int stabilizedMCS = 0; int numprocs = 1;

  double J = 1; double k = 1; double beta = 1/(k*T);
  double values[5]; long idum = -1; int count = 0; double mcsDouble, countDouble;

  ofstream file(filename, ofstream::binary);
  for (int mcs = mcsStart; mcs <= mcsStop; mcs+=mcsStep){
    solveGivenT(L, mcs, T, k, J, values, idum, stabilizedMCS, ordered, count,p);
    calculateVarNormalize(values, L, mcs, stabilizedMCS, numprocs);
    calculateCChi(values, k, T);
    mcsDouble = (double)mcs; countDouble = (double)count;
    file.write(reinterpret_cast<const char*>(&mcsDouble), sizeof(double));
    file.write(reinterpret_cast<const char*>(&T), sizeof(double));
    file.write(reinterpret_cast<const char*>(values), 5*sizeof(double));
    file.write(reinterpret_cast<const char*>(&countDouble), sizeof(double));
  }
  file.close();
}

void solveForDifferentMCS2(double T, bool ordered, string filename, int mcsStart, int mcsStop, int mcsStep){
  int L = 20; double p[1000]; int stabilizedMCS = 3e4; int numprocs = 1;

  double J = 1; double k = 1; double beta = 1/(k*T);
  double values[5]; long idum = -1; int count; double mcsDouble, countDouble;

  ofstream file(filename, ofstream::binary);
  for (int mcs = mcsStart; mcs <= mcsStop; mcs*=mcsStep){
    solveGivenT(L, mcs, T, k, J, values, idum, stabilizedMCS, ordered, count,p);
    calculateVarNormalize(values, L, mcs, stabilizedMCS, numprocs);
    calculateCChi(values, k, T);
    mcsDouble = (double)mcs; countDouble = (double)count;
    file.write(reinterpret_cast<const char*>(&mcsDouble), sizeof(double));
    file.write(reinterpret_cast<const char*>(&T), sizeof(double));
    file.write(reinterpret_cast<const char*>(values), 5*sizeof(double));
    file.write(reinterpret_cast<const char*>(&countDouble), sizeof(double));
  }
  file.close();
}

void calculateProbability(double T, string filename){
  int L = 20; int mcs = 1e6; bool ordered = false;
  double J = 1; double k = 1; double beta = 1/(k*T);
  double values[5]; long idum = -1; int count;
  double p[1000] = {0}; int stabilizedMCS = 3e4; int numprocs = 1;

  //ofstream file(filename, ofstream::binary);
  ofile.open(filename);

  solveGivenT(L, mcs, T, k, J, values, idum, stabilizedMCS, ordered, count,p);
  calculateVarNormalize(values, L, mcs, stabilizedMCS, numprocs);

  //file.write(reinterpret_cast<const char*>(&values[0]), sizeof(double));
  //file.write(reinterpret_cast<const char*>(&values[1]), sizeof(double));
  ofile << values[0] << ' ' << values[1] << endl;

  double energy;
  for (int i = 0; i<1000;++i){
    if (p[i] != 0){
      energy = -i/(double)(L*L);

      //file.write(reinterpret_cast<const char*>(&energy), sizeof(double));
      //file.write(reinterpret_cast<const char*>(&p[i]), sizeof(double));
      ofile << setw(5) << setprecision(3) << energy;
      ofile << setw(15) << setprecision(8) << p[i] << endl;
    }
  }
  ofile.close();
}
