/*
      PROGRAM FOR THE NUMERICAL SIMULATION
      OF A SCALAR FIELD THEORY IN D = 2
      USING THE PATH-INTEGRAL FORMULATION
      OF QUANTUM FIELD THEORY
*/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

// MACROS DEFINING THE LATTICE STRUCTURE
#define N_X 10
#define N_T 10
#define N_VOL (N_X * N_T)

// MACRO DEFINING PI GREGO
#define PIGR 3.141592654

// MACROS DEFINING THE RAN2 GENERATOR
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

typedef struct
{
      double field[N_X][N_T];
      int npp[N_X][2];
      int nmm[N_X][2];
} Lattice;

typedef struct
{
      int iflag;
      int measures;
      int i_decorrel;
      double extfield;
      double mass;
      double mass2;
      double mass2p4;
} Parameters;

typedef struct
{
      long idum;
      long idum2;
      long iv[NTAB];
      long iy;
      long seed;
} Ran2Generator;

void geometry(Lattice *pLattice);
void initialize_lattice(Lattice *pLattice, Parameters *pParameters, Ran2Generator *pRan2Generator);
void update_heatbath(Lattice *pLattice, Parameters *pParameters, Ran2Generator *pRan2Generator);
void update_overrelax(Lattice *pLattice, Parameters *pParameters);
void energy(Lattice *pLattice, Parameters *pParameters, FILE *pOut);
double ran2(Ran2Generator *pRan2Generator);
void ranstart(Ran2Generator *pRan2Generator);
void ranfinish(Ran2Generator *pRan2Generator);

int main()
{
      clock_t start = clock();

      Lattice lattice;
      Lattice *pLattice = NULL;
      pLattice = &lattice;

      Parameters parameters;
      Parameters *pParameters = NULL;
      pParameters = &parameters;

      Ran2Generator ran2Generator;
      Ran2Generator *pRan2Generator = NULL;
      pRan2Generator = &ran2Generator;

      FILE *pIn;
      FILE *pOut;
      FILE *pLatt;

      ranstart(pRan2Generator);
      
      // INPUT FILE CONTAINING THE PARAMETERS OF THE SIMULATION
      pIn = fopen("parameters.txt", "r");
      
      // OUTPUT FILE TO STORE THE MEASURES OF THE OBSERVABLES
      pOut = fopen("measures_out", "w"); 

      if(pIn == NULL) 
      {
            printf("Input file can't be opened.\n");
            return EXIT_FAILURE;
      }      

      if(fscanf(pIn, "%d %d %d %lf %lf", &parameters.iflag, 
            &parameters.measures, &parameters.i_decorrel, 
            &parameters.extfield, &parameters.mass) != 5)
      {
            printf("Error reading input file.\n");
      }

      parameters.mass2 = pow(parameters.mass, 2);
      parameters.mass2p4 = pow(parameters.mass, 2) + 4;
      fclose(pIn);

      // PRELIMINARY OPERATIONS
      geometry(pLattice);     // initialize boundary conditions
      initialize_lattice(pLattice, pParameters, pRan2Generator); // initial configuration
      // ----------------------------------------------------------

      for(int i = 0; i < parameters.measures; i++)
      {
            // UPDATING LATTICE CONFIGURATION ...
            // i_decorrel sweep of the whole lattice
            for(int j = 0; j < parameters.i_decorrel; j++)
            {
                  update_heatbath(pLattice, pParameters, pRan2Generator);
                  update_overrelax(pLattice, pParameters);
                  update_overrelax(pLattice, pParameters);
                  update_overrelax(pLattice, pParameters);
                  update_overrelax(pLattice, pParameters);
            }

            // MEASURES OF THE PHYSICAL OBSERVABLES
            energy(pLattice, pParameters, pOut);
      }

      fclose(pOut);

      // END OF MONTE-CARLO SIMULATION
      // ------------------------------------------------------------

      /*
      SAVING LATTICE CONFIGURATION AND RANDOM GENERATOR STATE
      TO POSSIBLY RESTART FROM THIS SITUATION
      */

      pLatt = fopen("lattice", "w");

      for(int i = 0; i < N_X; i++)
      {
            for(int j = 0; j < N_T; j++)
            {
                  fprintf(pLatt, "%lf\t", lattice.field[i][j]);
            }
            fprintf(pLatt, "\n");
      }

      fclose(pLatt);

      ranfinish(pRan2Generator);

      clock_t end = clock();

      double cpu_time = ((double) (end - start)) / CLOCKS_PER_SEC;

      printf("CPU execution time: %.0lf seconds", cpu_time);

      return EXIT_SUCCESS;
}

void geometry(Lattice *pLattice)
{
      /*
      npp and nmm are vectors of integers that take as input
      a coordinate and return the forward and backward
      coordinates, respectively, taking into account the 
      periodic boundary conditions.
      */

      /*
      The lattice may be anisotropic, i.e. the spatial and temporal
      direction may have different lengths (NX != NT)
      */

      for(int i = 0; i < N_X; i++)
      {
            pLattice->npp[i][1] = i + 1;
            pLattice->nmm[i][1] = i - 1;
      }
      pLattice->npp[N_X - 1][1] = 0;  // PERIODIC BOUNDARY CONDITIONS
      pLattice->nmm[0][1] = N_X - 1;  // AT THE EDGES OF THE LATTICE

      for(int i = 0; i < N_T; i++)
      {
            pLattice->npp[i][2] = i + 1;
            pLattice->nmm[i][2] = i - 1;
      }
      pLattice->npp[N_T - 1][2] = 0;  // PERIODIC BOUNDARY CONDITIONS
      pLattice->nmm[0][2] = N_T - 1;  // AT THE EDGES OF THE LATTICE
}

void initialize_lattice(Lattice *pLattice, Parameters *pParameters, Ran2Generator *pRan2Generator)
// SET UP THE STARTING CONFIGURATION OF THE MARKOV CHAIN
{
      FILE *pLatt;

      // COLD START ... (the field is set to zero as T = 0)
      if(pParameters->iflag == 0)
      {
            for(int i = 0; i < N_X; i++)
            {
                  for(int j = 0; j < N_T; j++)
                  {
                        pLattice->field[i][j] = 0.0;
                  }
            }
      }
      // ... HOT START ... (random field in [-1; 1] as T = infinite)
      else if(pParameters->iflag == 1)
      {
            for(int i = 0; i < N_X; i++)
            {
                  for(int j = 0; j < N_T; j++)
                  {
                        double x = ran2(pRan2Generator); // random number in [0; 1]
                        pLattice->field[i][j] = 1.0 - 2.0 * x;
                  }
            }
      }
      // ... OR STARTING AGAIN FROM THE PREVIOUS LATTICE CONFIGURATION
      else
      {
            pLatt = fopen("lattice", "r");
            if(pLatt == NULL)
            {
                  printf("Lattice file can't be found.\n");
                  return EXIT_FAILURE;
            }
            for(int i = 0; i < N_X; i++)
            {
                  for(int j = 0; j < N_T; j++)
                  {
                        if(fscanf(pLatt, "%lf", &(pLattice->field[i][j])) != 1)
                        {
                              printf("Error reading lattice file.\n");
                        }
                  }  
            }

            if(fscanf(pLatt, "%ld", &(pRan2Generator->seed)))
            {
                  printf("Error reading seed.\n");
            }

            fclose(pLatt);
      }      
}

void update_heatbath(Lattice *pLattice, Parameters *pParameters, Ran2Generator *pRan2Generator)
{
      for(int i = 0; i < N_X; i++)
      {
            for(int j = 0; j < N_T; j++)
            {
                  int ip = pLattice->npp[i][1];
                  int im = pLattice->nmm[i][1];
                  int jp = pLattice->npp[j][2];
                  int jm = pLattice->nmm[j][2];

                  double force = 0.0;
                  double phi = pLattice->field[i][j];

                  force = force + pLattice->field[ip][j];
                  force = force + pLattice->field[im][j];
                  force = force + pLattice->field[i][jp];
                  force = force + pLattice->field[i][jm];

                  double sigma2 = 1.0 / pParameters->mass2p4; // variance of the gaussian
                                                              // distribution is 1/(m^2 + 4)
                  double aver = force * sigma2; // average of the gaussian distribution
                                                // is force/(m^2 + 4)
                  
                  // BOX MULLER ALGORITHM
                  double x = log(ran2(pRan2Generator));
                  x = sqrt(sigma2) * sqrt(- 2.0 * x);
                  double y = x * cos(2.0 * PIGR * ran2(pRan2Generator)) + aver;
                  
                  pLattice->field[i][j] = y;
            }
      }

}

void update_overrelax(Lattice *pLattice, Parameters *pParameters)
{
      for(int i = 0; i < N_X; i++)
      {
            for(int j = 0; j < N_T; j++)
            {
                  int ip = pLattice->npp[i][1];
                  int im = pLattice->nmm[i][1];
                  int jp = pLattice->npp[j][2];
                  int jm = pLattice->nmm[j][2];

                  double force = 0.0;
                  double phi = pLattice->field[i][j];

                  force = force + pLattice->field[ip][j];
                  force = force + pLattice->field[im][j];
                  force = force + pLattice->field[i][jp];
                  force = force + pLattice->field[i][jm];

                  double aver = force / pParameters->mass2p4;
                  // average of the gaussian distribution is force/(m^2 + 4)

                  pLattice->field[i][j] = 2.0 * aver - phi;
            }
      }
}

void energy(Lattice *pLattice, Parameters *pParameters, FILE *pOut)
{
      double xene_mass = 0.0;
      double xene_spat = 0.0;
      double xene_temp = 0.0;
      double xene_tot = 0.0;

      for(int i = 0; i < N_X; i++)
      {
            for(int j = 0; j < N_T; j++)
            {
                  int ip = pLattice->npp[i][1];
                  int jp = pLattice->npp[j][2];

                  double phi = pLattice->field[i][j];
                  double force_s = pLattice->field[ip][j];
                  double force_t = pLattice->field[i][jp];

                  xene_mass = xene_mass + pParameters->mass2 * pow(phi, 2);
                  xene_spat = xene_spat - 2.0 * phi * force_s + 2.0 * pow(phi, 2);
                  xene_temp = xene_temp - 2.0 * phi * force_t + 2.0 * pow(phi, 2);
            }
      }

      xene_mass = xene_mass / ((double) N_VOL);
      xene_spat = xene_spat / ((double) N_VOL);
      xene_temp = xene_temp / ((double) N_VOL);
      xene_tot = xene_mass + xene_spat - xene_temp;

      fprintf(pOut, "%lf\t %lf\t %lf\t %lf", xene_mass, xene_spat, xene_temp, xene_tot);
}

/*
----------------------------------------------------------------------
      RANDOM NUMBER GENERATOR: standard ran2 from numerical recipes
----------------------------------------------------------------------
*/
double ran2(Ran2Generator *pRan2Generator)
/*
Long period (> 2 × 1018) random number generator of L’Ecuyer with Bays-Durham shuffle
and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 (exclusive of
the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter
idum between successive deviates in a sequence. RNMX should approximate the largest floating
value that is less than 1.
*/
//    ACHTUNG!! 
// La funzione originale trovata in Numerical Recipes prende come parametro
// d'ingresso long *idum. Nel codice, invece, vengono passati come parametri
// d'ingresso il puntatore dell'intera struct Ran2Generator, i cui valori
// idum, idum2, iv[NTAB], iy verranno opportunamente aggiornati nella funzione.
{
      int j;
      long k;
      // static long idum2=123456789;
      // static long iy=0;
      // static long iv[NTAB];
      double temp;
      if (pRan2Generator->idum <= 0) {
      if (-(pRan2Generator->idum) < 1) pRan2Generator->idum = 1;
      else pRan2Generator->idum = - (pRan2Generator->idum);
      pRan2Generator->idum2 = (pRan2Generator->idum);
      for (j = NTAB + 7; j >= 0; j--) {
      k = (pRan2Generator->idum) / IQ1;
      pRan2Generator->idum = IA1 * (pRan2Generator->idum - k * IQ1)- k * IR1;
      if (pRan2Generator->idum < 0) pRan2Generator->idum += IM1;
      if (j < NTAB) pRan2Generator->iv[j] = pRan2Generator->idum;
      }
      pRan2Generator->iy = pRan2Generator->iv[0];
      }
      k = (pRan2Generator->idum) / IQ1;
      pRan2Generator->idum = IA1 * (pRan2Generator->idum - k * IQ1) - k * IR1;
      if (pRan2Generator->idum < 0) pRan2Generator->idum += IM1;
      k = pRan2Generator->idum2 / IQ2;
      pRan2Generator->idum2 = IA2 * (pRan2Generator->idum2 - k * IQ2) - k * IR2;
      if (pRan2Generator->idum2 < 0) pRan2Generator->idum2 += IM2;
      j = pRan2Generator->iy / NDIV;
      pRan2Generator->iy = pRan2Generator->iv[j] - pRan2Generator->idum2;
      pRan2Generator->iv[j] = pRan2Generator->idum;
      if (pRan2Generator->iy < 1) pRan2Generator->iy += IMM1;
      /*
      Initialize.
      Be sure to prevent idum =0.
      Load the shuffle table (after 8 warm-ups).
      Start here when not initializing.
      Compute idum=(IA1*idum) % IM1 without
      overflows by Schrage’s method.
      Compute idum2=(IA2*idum) % IM2 likewise.
      Will be in the range 0..NTAB-1.
      Here idum is shuffled, idum and idum2 are
      combined to generate output.
      */
      if ((temp = AM * pRan2Generator->iy) > RNMX) return RNMX; //Because users don’t expect endpoint values.
      else return temp;
}

void ranstart(Ran2Generator *pRan2Generator)
{
      FILE *pSeed;
      pSeed = fopen("randomseed", "r");

      if(pSeed == NULL)
      {
            printf("Missing seed file.\n");
            return;
      }
      /*
      If the file "randomseed" doesn't exist, you have to create it and
      set the default value for idum, which is expected to be a negative number.
      */
      if(fscanf(pSeed, "%ld", &(pRan2Generator->idum)) != 1)
      {
            printf("Error reading idum. Use the default seed: idum = -1\n");
            fclose(pSeed);
            return;
      }

      if(fscanf(pSeed, "%ld", &(pRan2Generator->idum2)) != 1)
      {
            printf("No idum2 found. Initializating to default value.\n");
            pRan2Generator->idum2 = 123456789;
            fclose(pSeed);
      }

      for(int i = 0; i < NTAB; i++)
      {
           if(fscanf(pSeed, "%ld", &(pRan2Generator->iv[i])) != 1)
           {
                  printf("No iv[%d] found. Initializating to default value.\n", i);
                  pRan2Generator->iv[i] = 0;
                  fclose(pSeed);
           } 
      }

      if(fscanf(pSeed, "%ld", &(pRan2Generator->iy)) != 1)
      {
            printf("No iy found. Initializating to default value.\n");
            pRan2Generator->iy = 0;
            fclose(pSeed);
      }

      if(pRan2Generator->idum >= 0)
      {
            pRan2Generator->idum = - pRan2Generator->idum - 1;
      }

      fclose(pSeed);
}

void ranfinish(Ran2Generator *pRan2Generator)
{
      FILE *pSeed;
      pSeed = fopen("randomseed", "w");

      fprintf(pSeed, "%ld\n", pRan2Generator->idum);

      fprintf(pSeed, "%ld\n", pRan2Generator->idum2);

      for(int i = 0; i < NTAB; i++)
      {
            fprintf(pSeed, "%ld\n", pRan2Generator->iv[i]);
      }
      
      fprintf(pSeed, "%ld\n", pRan2Generator->iy);

      fclose(pSeed);
}