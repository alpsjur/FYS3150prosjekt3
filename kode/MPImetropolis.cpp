/*
Parallellisert Ising-modell.

Kjør med 2 threads:
Project-path username$ mpirun -n 2 ./MPImetropolis.exe filnavn
*/

#include "metropolis.hpp"
#include "mpi.h"

ofstream ofile;
void writeToFile(double *values, double T);
void solveGivenT(int L, int mcs, double T, double k, double J, double *values,
                 long &idum, bool ordered);

int main(int argc, char* argv[]){
  int L = 2;
  int mcs = 1e5;
  double J = 1;
  double k = 1;
  double initialT = 1;
  double finalT = 2;
  double dT = 0.1;
  bool ordered = true;


  char *outfilename;
  long idum;
  int my_rank, numprocs;
  double values[5], collectedValues[5], E, M;
  double p[1000];

  //  MPI initsialisering
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

  if (my_rank == 0) {
    outfilename = argv[1];
    ofile.open(outfilename);
  }


  //bestemmer hvor mange mcs hver prosess skal kjøre.
  int my_mcs = mcs/numprocs;
  if ( (my_rank == numprocs-1) &&((mcs%numprocs)!=0) ) my_mcs += mcs%numprocs;

  //sender felles variabler til alle noder
  MPI_Bcast (&L, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast (&initialT, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&finalT, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast (&dT, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  //hver node har et eget tilfeldig-tall-frø/seed
  idum = -1-my_rank;  // random starting point

  //looper over alle T-verdien, tar tiden
  double  TimeStart, TimeEnd, TotalTime;
  TimeStart = MPI_Wtime();
  for ( double T = initialT; T <= finalT; T += dT){
    int count = 0;
    solveGivenT(L, my_mcs, T, k, J, values, idum, ordered);

    // finner totalt gjennomsnitt
    for( int i =0; i < 5; i++){
      MPI_Reduce(&values[i], &collectedValues[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    // printer resultat
    if ( my_rank == 0) {
      calculateVarNormalize(collectedValues, L, mcs);
      calculateCChi(collectedValues, k, T);
      writeToFile(collectedValues, T);
      //output(L, mcs, T, collectedValues);
    }
  }


  ofile.close();  // close output file
  TimeEnd = MPI_Wtime();
  TotalTime = TimeEnd-TimeStart;
  if ( my_rank == 0) {
    cout << "Time = " <<  TotalTime  << " on number of processors: "  << numprocs  << endl;
  }

  // End MPI
  MPI_Finalize ();

  return 0;
}

void solveGivenT(int L, int mcs, double T, double k, double J, double *values,
                 long &idum, bool ordered){
  double beta = 1/(k*T);
  double w[17];
  imat spinMatrix;
  double nSpins = L*L;
  int sumP = 0;
  int deltaE;
  //regner ut mulige w-verdier
  for (int i=0; i<5; ++i){
    w[i*4] = exp(-beta*J*(i*4-8));
  }
  double E = 0; double M = 0;
  initialize(L, spinMatrix, values, E, M, J, ordered);
  //går gjennom gitt antall monte carlo sykluser (mcs)
  for (int i = 0; i < mcs; ++i){
    //går gjennom alle spinnene
    for (int j=0; j < nSpins; ++j){
      //velger et tilfeldig spinn
      int k = (int) (ran2(&idum)*(double)L);
      int l = (int) (ran2(&idum)*(double)L);
      //beregner endring i energi
      deltaE = (spinMatrix(k,periodic(l+1,L))+
                spinMatrix(k,periodic(l-1,L))+
                spinMatrix(periodic(k+1,L),l)+
                spinMatrix(periodic(k-1,L),l))*
                spinMatrix(k,l)*2;
      //Utfører metropolis-testen
      if (ran2(&idum) <= w[deltaE+8]){
        spinMatrix(k,l) *= -1;
        E += (double) deltaE*J;
        M += (double) 2*spinMatrix(k,l);
      }
    }
    //metropolis(L, spinMatrix, idum, E, M, w, J, count, p, i);
    values[0] += E;
    values[1] += E*E;
    values[2] += fabs(M);
    values[3] += M*M;
    values[4] += M;
  }
}


void writeToFile(double *values, double T){
  ofile << setiosflags(ios::showpoint | ios::uppercase);
  ofile << setw(15) << setprecision(8) << T;
  ofile << setw(15) << setprecision(8) << values[0];
  ofile << setw(15) << setprecision(8) << values[1];
  ofile << setw(15) << setprecision(8) << values[2];
  ofile << setw(15) << setprecision(8) << values[3];
  ofile << setw(15) << setprecision(8) << values[4] << endl;
}
