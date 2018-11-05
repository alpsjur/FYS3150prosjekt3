#include "catch.hpp"
#include "metropolis.hpp"

TEST_CASE("Tester om vi faar riktige forventningsverdier for L=2"){
  double J = 1; double T = 1; double k = 1; int L = 2;
  double beta = 1/(k*T); int mcs = 1e7;
  double values[5]; long idum = -1;
  int nSpins = L*L; bool ordered=true;

  solveGivenT(L, mcs, T, k, J, values, idum, ordered);
  calculateVarNormalize(values, L, mcs);
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
  analyticalE /= Z; analyticalM /= Z;
  EE /= Z; MM = MM /= Z;
  double analyticalC = (EE-analyticalE*analyticalE)/(k*T*T*nSpins);
  double analyticalChi = (MM-analyticalM*analyticalM)/(k*T*nSpins);
  analyticalE /= nSpins;
  analyticalM /= nSpins;
  analyticalAbsM /= nSpins*Z;

  //cout << analyticalE << " " << analyticalC << " " << analyticalM << " " << analyticalChi << endl;
  REQUIRE(fabs(values[0]/analyticalE) == Approx(1).epsilon(0.001));
  REQUIRE(fabs(values[1]/analyticalC) == Approx(1).epsilon(0.01));
  REQUIRE(fabs(values[2]/analyticalAbsM) == Approx(1).epsilon(0.001));
  REQUIRE(fabs(values[3]/analyticalChi) == Approx(1).epsilon(0.01));

}
