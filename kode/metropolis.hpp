#ifndef METROPOLIS_H
#define	METROPOLIS_H

//#include "catch.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <armadillo>
#include <string>
using namespace  std;
using namespace arma;

//inline function som bestemmer periodiske grenser
inline int periodic(int i, int L){
  return (i+L)%(L);
}

//funksjon som beregner total energi for systemet
double calculateTotalEnergy(int, imat, double);

//funksjon som initsialiserer systemet
void initialize(int L, imat &spinMatrix, double *values,double &E, double &M, double J, bool);

//finner forventningsverdier gitt temperatur ved MonteCarlo-metoden
void solveGivenT(int, int, double, double, double, double *, long &, bool, int&, double *);

//bergener variansen til E og M, og nrmaliserer med hensyn på antall monte carlo-sykluser og antall spinn
void calculateVarNormalize(double *values, int L, int mcs);

//funnksjon som beregner varmekapasiteten og susceptibiliteten
void calculateCChi(double *values, double k, double T);

//skriver beregnede verdier til fil
void writeToFile(ofstream ofile,double *values, double T, int mcs);

//parallelisert løsning for flere T-verdier
void MPIsolve(int argc, char* argv[], string outfilename, int L, int mcs,
              double J, double k ,double initialT, double finalT, double dT,
              bool ordered);

//random number generator, initsialiseres med negativt frø/seed. Kopiert fra Mortens kode
double ran2(long*);

#endif /* METROPOLIS_H */
