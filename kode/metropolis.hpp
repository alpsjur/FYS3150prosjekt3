#ifndef METROPOLIS_H
#define	METROPOLIS_H

//#include "catch.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <armadillo>
using namespace  std;
using namespace arma;

//inline function som bestemmer periodiske grenser
inline int periodic(int i, int L){
  return (i+L)%(L);
}

//funksjon som beregner total energi for systemet
double calculateTotalEnergy(int, imat, double);

//funksjon som initsialiserer systemet
void initialize(int L, imat &spinMatrix, double &E, double &M, double J, bool);

//mertropolis-algoritmen
void metropolis(int L, imat &spinMatrix, long &idum, double &E, double &M,double *w, double J);

//
void solveGivenT(int, int, double, double, double, double *, long &, bool);

//random number generator, initsialiseres med negativt frø/seed
double ran2(long*);

#endif /* METROPOLIS_H */
