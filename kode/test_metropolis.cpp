#include "catch.hpp"
#include "metropolis.hpp"

TEST_CASE("Tester om vi faar riktige forventningsverdier når L=2"){
  double J = 1; double T = 1; double k = 1; int L = 2;
  double beta = 1/(k*T); int mcs = 1e7;
  double values[5]; long idum = -1;
  int nSpins = L*L; bool ordered=true; int count = 0;
  double numP[1500]; int stabilizedMCS = 1e4; int numprocs = 1;

  solveGivenT(L, mcs, T, k, J, values, idum, stabilizedMCS, ordered, count, numP);
  calculateVarNormalize(values, L, mcs, stabilizedMCS, numprocs);
  calculateCChi(values, k, T);


  int degeneracy[6] = {1,4,4,2,4,1};
  int energy[6] = {-8,0,0,8,0,-8};
  int magnetization[6] = {4,2,0,0,-2,-4};

  double p;
  double Z = 0;
  double analyticalE = 0; double analyticalAbsM = 0; double analyticalM = 0;
  double EE = 0;  double MM = 0;
  for (int i=0; i<6; ++i){
    p = degeneracy[i]*exp(-beta*energy[i]);
    Z += p;
    analyticalE += energy[i]*p;
    analyticalM += magnetization[i]*p;
    analyticalAbsM += fabs(magnetization[i])*p;
    EE += energy[i]*energy[i]*p;
    MM += magnetization[i]*magnetization[i]*p;
  }
  analyticalE /= Z; analyticalM /= Z; analyticalAbsM /= Z;
  EE /= Z; MM = MM /= Z;
  double analyticalC = (EE-analyticalE*analyticalE)/(k*T*T*nSpins);
  double analyticalChi = (MM-analyticalAbsM*analyticalAbsM)/(k*T*nSpins);
  analyticalE /= nSpins;
  analyticalM /= nSpins;
  analyticalAbsM /= nSpins;

  //cout << values[0] << values [1] << values[2] << values[3] << endl;
  //cout << analyticalE << " " << analyticalC << " " << analyticalAbsM << " " << analyticalChi << endl;
  REQUIRE(fabs(values[0]/analyticalE) == Approx(1).epsilon(0.001));
  REQUIRE(fabs(values[1]/analyticalC) == Approx(1).epsilon(0.01));
  REQUIRE(fabs(values[2]/analyticalAbsM) == Approx(1).epsilon(0.001));
  REQUIRE(fabs(values[3]/analyticalChi) == Approx(1).epsilon(0.01));

}

TEST_CASE("Tester om summen av beregned sannsynlighet er 1"){
  double J = 1; double k = 1; int L = 10;
  int mcs = 1e6;
  double values[5]; long idum = -1;
  int nSpins = L*L; bool ordered=true; int count = 0;
  double probability1[1000] = {0}; double probability2[1000] = {0};
  double sumProbability1 = 0; double sumProbability2 = 0;
  int stabilizedMCS = 3e4;

  solveGivenT(L, mcs, 1, k, J, values, idum, stabilizedMCS, ordered, count, probability1);
  solveGivenT(L, mcs, 2.4, k, J, values, idum, stabilizedMCS, ordered, count, probability2);
  for (int i = 0; i<1000;++i){
    sumProbability1 += probability1[i];
    sumProbability2 += probability2[i];
  }
  REQUIRE(sumProbability1 == Approx(1));
  REQUIRE(sumProbability2 == Approx(1));
}
