//#include "mpi.h"
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
  return (i+L)%L;
}

//funksjon som beregner total energi for systemet
double calculateTotalEnergy(int, imat);

//funksjon som initsialiserer systemet
void initialize(int, imat&, double&);


int main(int argc, char const *argv[]) {

  return 0;
}

void MonteCarlo(int nSpins){
  //itererer over alle spinnene
  for (int i=0; i < nSpins; ++k){

  }
  return;
}

void initialize(int L, imat &spinMatrix, double &E){
  //lager uniform initsialtilstand
  //spinMatrix.ones();

  //lager tilfeldig initsialtilstand
  spinMatrix = randi<imat>(L, L, distr_param(0, 1));   //setter matriseelementene til 0 eller 1
  for (int k = 0; k < L; ++k){
    for (int l = 0; l< L; ++l){
      if (spinMatrix(k,l) == 0){                       //elementene som er 0 endres til -1
        spinMatrix(k,l) = -1;
      }
    }
  }
  E = calculateTotalEnergy(L, spinMatrix);
  return;
}

double calculateTotalEnergy(int L, imat spinMatrix){
  //kalkulerer total energi for hele systemet ved å iterere over alle elementene
  double E = 0;
  for (int k=0; k < L; ++k){
    for (int l=0; l < L; ++l){
      E -= (spinMatrix(k,periodic(l-1,L))+
            spinMatrix(periodic(k-1,l),l))*
            spinMatrix(k,l);
    }
  }
  return E;
}
