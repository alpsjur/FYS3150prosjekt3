/*
Parallellisert Ising-modell.

Kjør med 2 threads:
Project-path username$ mpirun -n 2 ./MPImetropolis.exe filnavn
*/

#include "metropolis.hpp"
#include "mpi.h"


//ofstream ofile;
void writeToFile(double *values, double T);
void solveGivenT(int L, int mcs, double T, double k, double J, double *values,
                 long &idum, int, bool ordered);

int main(int argc, char* argv[]){
  int L = atoi(argv[2]);
  int mcs = 1e6;
  double J = 1;
  double k = 1;
  double initialT = 2.20;
  double finalT = 2.30;
  double dT = 0.005;
  bool ordered = false;
  int stabilizedMCS = 25e3;


  char *outfilename;
  long idum;
  int my_rank, numprocs;
  double values[5], collectedValues[5];

  //  MPI initsialisering
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

  ofstream ofile;
  if (my_rank == 0) {
    outfilename = argv[1];
    ofile.open(outfilename, ofstream::binary);
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
    solveGivenT(L, my_mcs, T, k, J, values, idum, stabilizedMCS, ordered);

    // finner totalt gjennomsnitt
    for( int i =0; i < 5; i++){
      MPI_Reduce(&values[i], &collectedValues[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    // printer resultat
    if ( my_rank == 0) {
      calculateVarNormalize(collectedValues, L, mcs, stabilizedMCS,numprocs);
      calculateCChi(collectedValues, k, T);
      ofile.write(reinterpret_cast<const char*>(&T), sizeof(double));
      ofile.write(reinterpret_cast<const char*>(collectedValues), 5*sizeof(double));
      cout << "Ferdig med temperatur " << T << endl;
    }
  }


  TimeEnd = MPI_Wtime();
  TotalTime = TimeEnd-TimeStart;
  if ( my_rank == 0) {
    ofile.close();  // close output file
    cout << "Time = " <<  TotalTime  << " on number of processors: "  << numprocs  << endl;
  }

  // End MPI
  MPI_Finalize ();
  return 0;
}

void solveGivenT(int L, int mcs, double T, double k, double J, double *values,
                 long &idum, int stabilizedMCS, bool ordered){
  double beta = 1/(k*T);
  double w[17];
  imat spinMatrix;
  double nSpins = L*L;
  int deltaE, sumP = 0;
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
      int l = (int) (ran2(&idum)*(double)L);
      int m = (int) (ran2(&idum)*(double)L);
      //beregner endring i energi
      deltaE = (spinMatrix(l,periodic(m+1,L))+
                spinMatrix(l,periodic(m-1,L))+
                spinMatrix(periodic(l+1,L),m)+
                spinMatrix(periodic(l-1,L),m))*
                spinMatrix(l,m)*2;
      //Utfører metropolis-testen
      if (ran2(&idum) <= w[deltaE+8]){
        spinMatrix(l,m) *= -1;
        E += (double) deltaE*J;
        M += (double) 2*spinMatrix(l,m);
      }
    }
    if (i > stabilizedMCS){
      values[0] += (double) E;
      values[1] += (double) E*E;
      values[2] += (double) fabs(M);
      values[3] += (double) M*M;
      values[4] += (double) M;
    }
  }
}
