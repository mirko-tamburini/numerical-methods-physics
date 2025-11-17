/*
      PROGRAM FOR THE NUMERICAL SIMULATION
      OF THE HARMONIC OSCILLATOR IN D = 1
      USING THE PATH-INTEGRAL FORMULATION
      OF QUANTUM MECHANICS
*/

#include <stdio.h>
#include <math.h>
#include <time.h>

// MACROS DEFINING THE LATTICE STRUCTURE
#define N_LATT 100

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
      double field[N_LATT];
      int npp[N_LATT];
      int nmm[N_LATT];
} Lattice;

typedef struct
{
      int iflag;
      int measures;
      int i_decorrel;
      double i_term;
      double eta;
      double d_metro;
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
void measure(Lattice *pLattice, FILE *pOut);
void update_metropolis(Lattice *pLattice, Parameters *pParameters, Ran2Generator *pRan2Generator);
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
            return 0;
      }      

      if(fscanf(pIn, "%d %d %d %lf %lf", &parameters.iflag, 
            &parameters.measures, &parameters.i_decorrel, 
            &parameters.i_term, &parameters.eta) != 5)
      {
            printf("Error reading input file.\n");
      }

      // if(fscanf(pIn, "%lf", &parameters.d_metro) !=1)
      // {
      //       printf("Error reading d_metro in input file.\n");
      // }

      parameters.d_metro = 2.0 * sqrt(parameters.eta);

      fclose(pIn);

      // PRELIMINARY OPERATIONS
      geometry(pLattice);     // initialize boundary conditions
      initialize_lattice(pLattice, pParameters, pRan2Generator); // initial configuration
      // ----------------------------------------------------------

      // THERMALIZATION
      for(int i = 0; i < parameters.i_term; i++)
      {
            update_metropolis(pLattice, pParameters, pRan2Generator);
      }
      // -----------------------------------------------------------

      // IN-EQUILIBRIUM SESSION WITH MEASURES
      for(int i = 0; i < parameters.measures; i++)
      {
            // UPDATING LATTICE CONFIGURATION ...
            for(int j = 0; j < parameters.i_decorrel; j++)
            {
                  update_metropolis(pLattice, pParameters, pRan2Generator);
            }

            measure(pLattice, pOut);

      }

      fclose(pOut);

      // END OF MONTE-CARLO SIMULATION

      /*
      SAVING LATTICE CONFIGURATION AND RANDOM GENERATOR STATE
      TO POSSIBLY RESTART FROM THIS SITUATION
      */

      pLatt = fopen("lattice", "w");

      for(int i = 0; i < N_LATT; i++)
      {
            fprintf(pLatt, "%lf\n", lattice.field[i]);
      }

      fclose(pLatt);

      ranfinish(pRan2Generator);

      clock_t end = clock();

      double cpu_time = ((double) (end - start)) / CLOCKS_PER_SEC;

      printf("CPU execution time: %.0lf seconds", cpu_time);

      return 0;
}

void geometry(Lattice *pLattice)
{
      /*
      npp and nmm are vectors of integers that take as input
      a coordinate and return the forward and backward
      coordinates, respectively, taking into account the 
      periodic boundary conditions.
      */

      for(int i = 0; i < N_LATT; i++)
      {
            pLattice->npp[i] = i + 1;
            pLattice->nmm[i] = i - 1;
      }
      pLattice->npp[N_LATT - 1] = 0;  // PERIODIC BOUNDARY CONDITIONS
      pLattice->nmm[0] = N_LATT - 1;  // AT THE EDGES OF THE LATTICE
}

void initialize_lattice(Lattice *pLattice, Parameters *pParameters, Ran2Generator *pRan2Generator)
{
      FILE *pLatt;

      // COLD START ... (T = 0)
      if(pParameters->iflag == 0)
      {
            for(int i = 0; i < N_LATT; i++)
            {
                  pLattice->field[i] = 0.0;
            }
      }
      // ... HOT START ... (T = infinite)
      else if(pParameters->iflag == 1)
      {
            for(int i = 0; i < N_LATT; i++)
            {
                  double x = ran2(pRan2Generator);
                  x = 1.0 - 2.0 * x;      // random number in [-1; 1]
                  pLattice->field[i] = x; 
            }
      }
      // ... OR STARTING AGAIN FROM THE PREVIOUS LATTICE CONFIGURATION
      else
      {
            pLatt = fopen("lattice", "r");
            if(pLatt == NULL)
            {
                  printf("Lattice file can't be found.\n");
            }
            for(int i = 0; i < N_LATT; i++)
            {
                  if(fscanf(pLatt, "%lf", &(pLattice->field[i])) != 1)
                  {
                        printf("Error reading lattice file.\n");
                  }
            }
            // if(fscanf(pLatt, "%ld", &(pRan2Generator->seed)))
            // {
            //       printf("Error reading seed.\n");
            // }
            fclose(pLatt);
      }      
}

void update_metropolis(Lattice *pLattice, Parameters *pParameters, Ran2Generator *pRan2Generator)
{
      double c1 = 1.0 / pParameters->eta;
      double c2 = (1.0 / pParameters->eta + pParameters->eta / 2.0);

      for(int i = 0; i < N_LATT; i++)
      {
            int ip = pLattice->npp[i];
            int im = pLattice->nmm[i];

            double force = pLattice->field[ip] + pLattice->field[im];

            double phi = pLattice->field[i];
            double phi_test = phi + 2.0 * pParameters->d_metro * (0.5 - ran2(pRan2Generator));
            
            double p_rat = c1 * phi_test * force - c2 * pow(phi_test, 2);
            p_rat = p_rat - c1 * phi * force + c2 * pow(phi, 2);
            
            double x = ran2(pRan2Generator);
            x = log(x);
            
            if(x < p_rat)
            {
                  pLattice->field[i] = phi_test;
            }
      }
}

void measure(Lattice *pLattice, FILE *pOut)
{
      double obs1 = 0.0;
      double obs2 = 0.0;

      for(int i = 0; i < N_LATT; i++)
      {
            int ip = pLattice->npp[i];
            obs1 = obs1 + pow(pLattice->field[i], 2);
            obs2 = obs2 + pow((pLattice->field[i] - pLattice->field[ip]), 2);
      }

      obs1 = obs1/N_LATT;
      obs2 = obs2/N_LATT;

      fprintf(pOut, "%lf %lf\n", obs1, obs2);
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