/***************************************************************************
 *   Copyright (C) 2007 by vahid mokhtari and Ramin Fathzadeh              *
 *   mokhtari@mrl.ir                                                       *
 *   fathzadeh@mrl.ir                                                      *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include "efcm.h"

/***************************************************************************/
/*                            Fuzzy C-Means Class                          */
/***************************************************************************/

/* Constructors ===========================================================*/

FuzzyCMeans::FuzzyCMeans()
{
        Epsilon                 = 0.00001;
        m                       = 2.0;
        C                       = 2;
        D                       = 2;
        N                       = 2;
        NumberOfIterations      = 0;
        K                       = 2;
        S                       = 0;
        Validation              = 0;
        X                       = NULL;
        V                       = NULL;
        U                       = NULL;
        MaxValue                = NULL;
        ConfusionMatrix         = NULL;
}

FuzzyCMeans::FuzzyCMeans( double epsilon, double fuzzyexp, int c, int d,
                          long n, int tc, int v, int s, double *maxvalue )
{
        Epsilon                 = epsilon;
        m                       = fuzzyexp;
        C                       = c;
        D                       = d;
        N                       = n;
        NumberOfIterations      = 0;
        K                       = tc;
        S                       = s;
        Validation              = v;
        X                       = NULL;
        V                       = NULL;
        U                       = NULL;
        MaxValue                = maxvalue;
        ConfusionMatrix         = NULL;
}

FuzzyCMeans::~ FuzzyCMeans()
{
        if( V != NULL )
        {
                for( int i = 0; i < C; i++ )
                        delete [] V[i];
                delete [] V;
        }
        if( U != NULL )
        {
                for( long i = 0; i < N; i++ )
                        delete [] U[i];
                delete [] U;
                U = NULL;
        }
}

/* Parameters =============================================================*/

void FuzzyCMeans::ResetParameters( double epsilon, double fuzzyexp, int c, int d, long n, double *maxvalue )
{
        Epsilon                 = epsilon;
        m                       = fuzzyexp;
        C                       = c;
        D                       = d;
        N                       = n;
        NumberOfIterations      = 0;
        X                       = NULL;
        V                       = NULL;
        U                       = NULL;
        MaxValue                = maxvalue;
}

void FuzzyCMeans::GetParameters( double &epsilon, double &fuzzyexp, int &c, int &d, long &n )
{
        epsilon                 = Epsilon;
        fuzzyexp                = m;
        c                       = C;
        d                       = D;
        n                       = N;
}

/* Main Methods ===========================================================*/

int FuzzyCMeans::FCM( int displayinfo )
{
    double SqrError = 2 * Epsilon;
   
    /* Dynamically create various data structures */
    Init();

    /* Run the updates iteratively */
    while( SqrError > Epsilon )
    {
           NumberOfIterations++;
           UpdateCentroids();
           SqrError = UpdateUmatrix();
           if( displayinfo == 1 )
               printf("\n\titeration(%d) square error:(%f)",NumberOfIterations, SqrError);
    }

    /* We go ahead and update the centroids - presumably this will not 
       change much, since the overall square error in U is small */
    UpdateCentroids();

   return 0;
}

/*=========================================================================*/

/* Allocate storage for U and V dynamically. */
int FuzzyCMeans::Init()
{
  if( X == NULL )
  {
      fprintf(stderr,"%s:%d\nNo storage allocated to **X... [Fail]\n", __FILE__, __LINE__);
      exit(1);
  }
        
  long i;
  int  j;

  NumberOfIterations = 0;

  /* Allocate necessary storage */
  V = new double * [C];
  for( j = 0; j < C; j++ )
    V[j] = new double [D];

  U = new double * [N];
  for( i = 0; i < N; i++ )
    U[i] = new double [C];

  /* Place random values in V, then update U matrix based on it */
  if( MaxValue != NULL )
  {
      if( OS == 0 )
      {
        srand((unsigned) time(NULL));
        for( i = 0; i < C; i++ )
            for( j = 0; j < D; j++ )
                V[i][j] = (rand() * 0.00001) * MaxValue[j];
      }
      else
      {
        #ifdef _UNIX_
        srand48((unsigned) time(NULL));
        for( i = 0; i < C; i++ )
            for( j = 0; j < D; j++ )
                V[i][j] = drand48() * MaxValue[j];
        #endif
      }
  }
  
  /* Once values are populated in V, update the U matrix to save values */
  UpdateUmatrix();

  return 0;
}

/*=========================================================================*/

/* 
   UpdateCentroids()
    Given a membership matrix U, recalculate the cluster centroids as the
    "weighted" mean of each contributing example from the dataset. Each
    example contributes by an amount proportional to the membership value.
*/
int FuzzyCMeans::UpdateCentroids()
{
  long   k;
  int    i, j;
  double Numerator[D], Denominator;
  double U_ikm;

  /* For each cluster */
  for( i = 0; i < C; i++ )
  {
    /* Zero out numerator and denominator options */
    Denominator = 0;
    for( j = 0; j < D; j++ ) 
      Numerator[j] = 0;

    /* Calculate numerator and denominator together */
    for( k = 0; k < N; k++ )
    {
      U_ikm = pow( U(i,k) , m );
      Denominator += U_ikm;
      for( j = 0; j < D; j++ ) 
        Numerator[j] += U_ikm * X[k][j];
    }

    /* Calculate V */
    for( j = 0; j < D; j++ )
      V[i][j] = Numerator[j] / Denominator;

  }  /* endfor: C clusters */

  return 0;
}

/*=========================================================================*/

double FuzzyCMeans::UpdateUmatrix()
{
  long   k;
  int    i, j;
  int    ExampleIsCentroid;
  double Summation, D_k[C];
  double SquareDifference = 0;
  double NewU;

  /* For each example in the dataset */
  for( k = 0; k < N; k++ )
  {
    /* Special case: If Example is equal to a Cluster Centroid,
       then U=1.0 for that cluster and 0 for all others */
    if( (ExampleIsCentroid = IsExampleCentroid(k)) != -1 )
    {
      for( i = 0; i < C; i++ )
      {
        if( i == ExampleIsCentroid )
          U(i,k) = 1.0;
        else 
          U(i,k) = 0.0;
      }
      continue;
    }

    /* Cache the distance between this vector and all centroids. */
    for( i = 0; i < C; i++ ) 
      D_k[i] = Distance( X[k] , V[i] );
    
    /* For each class */
    for( i = 0; i < C; i++ )
    {
      Summation = 0;

      /* Calculate summation */
      for( j = 0; j < C; j++ )
      {
        if( i == j ) 
          Summation += 1.0;
        else
          Summation += pow( D_k[i] / D_k[j] , (2.0/ (m-1)) );
      }

      /* Weight is 1/sum */
      NewU = 1.0 / (double)Summation;
      
      /* Add to the squareDifference */
      SquareDifference += (U(i,k) - NewU) * (U(i,k) - NewU);
      
      U(i,k) = NewU;
    }

  } /* endfor n */

  return SquareDifference;
}

/*=========================================================================*/

double FuzzyCMeans::ClusterValidation( )
{
    long   i;
    int    j;
    int    F[N];

    /* Allocate necessary storage to **Clusters */
    ConfusionMatrix = new long * [K];
    for( i = 0; i < K; i++ ) 
       ConfusionMatrix[i] = new long [C];

    for( i = 0; i < K; i++ )
      for( j = 0; j < C; j++ ) 
         ConfusionMatrix[i][j] = 0;

    /* set final cluster for every data */
    int max;
    for( i = 0; i < N; i++ )
    {
      for( max = j = 0; j < C; j++ )
           if ( U[i][j] > U[i][max] ) max = j;
      F[i] = max;
    }

    for( i = 0; i < K; i++ )
      for( j = 0; j < C; j++ )
      {
           for( long x = 0; x < N; x++ )
             if( TC[x] == i + S )
             {
               for( long y = 0; y < N; y++ )
                 if( F[y] == j )
                   if( x == y )
                     ConfusionMatrix[i][j]++;
             }
      }

    long double N11, N10, N01, N00;
    /* The category N11 contains the pairs of points that are in the same 
       cluster both in C and in C'. The category N10 contains the pairs of 
       points that are in the same cluster in C but not in C'. The 
       definitions of N01 and N00 are similar. */

    N11 = N10 = N01 = N00 = 0;

    // ----------N11
    for( i = 0; i < K; i++ )
      for( j = 0; j < C; j++ )
        N11 += ConfusionMatrix[i][j] * ConfusionMatrix[i][j];
    N11 -= N;
    N11 /= 2;

    // ----------N10
    long n = 0;
    for( j = 0; j < C; j++ )
    {
      n = 0;
      for( i = 0; i < K; i++ )
         n += ConfusionMatrix[i][j];
      N10 += n * n;
    }

    n = 0;
    for( i = 0; i < K; i++ )
      for( j = 0; j < C; j++ )
         n += ConfusionMatrix[i][j] * ConfusionMatrix[i][j];

    N10 -= n;
    N10 /= 2;

    // ----------N01
    n = 0;
    for( i = 0; i < K; i++ )
    {
      n = 0;
      for( j = 0; j < C; j++ )
         n += ConfusionMatrix[i][j];
      N01 += n * n;
    }

    n = 0;
    for( i = 0; i < K; i++ )
      for( j = 0; j < C; j++ )
         n += ConfusionMatrix[i][j] * ConfusionMatrix[i][j];

    N01 -= n;
    N01 /= 2;

    // ----------N00
    long double PP;
    double Jaccard;
    /* PP = N(N - 1)/2 is the total number of point pairs */

    PP = N * ( N - 1 ) / 2;
    N00 = PP - N11 - N10 - N01;
    
    Jaccard = N11 / ( N11 + N01 + N10 );
    
    for( i = 0; i < K; i++ )
         delete [] ConfusionMatrix[i];
    delete [] ConfusionMatrix;
    ConfusionMatrix = 0;

    return Jaccard;
}

/* Utilities ==============================================================*/

/* If X[k] == V[i] for some i, then return that i. Otherwise, return -1 */
int FuzzyCMeans::IsExampleCentroid(int k)
{
  int  i, j;
  
  for( i = 0; i < C; i++ )
  {
    for( j = 0; j < D; j++ )
    {
      if( X[k][j] != V[i][j] ) break;
    }
    if( j == D )  /* X==V */
      return i;
  }

  return -1;
}

/*=========================================================================*/

double FuzzyCMeans::Distance(double *v1, double *v2)
{
  int    x;
  double Sum = 0;
  
  for( x = 0; x < D; x++ ) 
    Sum += (v1[x]-v2[x]) * (v1[x]-v2[x]);
    
  return sqrt(Sum);
}

/* Loading Data ===========================================================*/

void FuzzyCMeans::SetData( double **x )
{
        X = x;
}

/*=========================================================================*/

void FuzzyCMeans::SetTrueClusters( int * tc )
{
        TC = tc;
}

/*=========================================================================*/

void FuzzyCMeans::SetMaxValue( double * maxvalue )
{
        MaxValue = maxvalue;
}

/* Output Utilities =======================================================*/

double ** FuzzyCMeans::GetCentroids()
{
  return V;
}

/*=========================================================================*/

double ** FuzzyCMeans::GetUmatrix()
{
  return U;
}

/*=========================================================================*/

int * FuzzyCMeans::GetMembers( int * M )
{
  long i;
  int  j, max;

  for( i = 0; i < N; i++ )
  {
    for( max = j = 0; j < C; j++ )
      if ( U[i][j] > U[i][max] ) max = j;
    M[i] = max;
  }

  return M;
}

/*=========================================================================*/

int FuzzyCMeans::GetNrFinalCluster()
{
    return C;
}

/*=========================================================================*/

int FuzzyCMeans::GetFinalClusters( int * f )
{
    int max, j;
    for( long i = 0; i < N ;i++ )
    {
      for( max = j = 0; j < C; j++ )
        if ( U[i][j] > U[i][max] ) max = j;
      *(f+i) = max;
    }
    return C;
}

/*=========================================================================*/

int FuzzyCMeans::GetClusters( double **v )
{
    for( int i = 0; i < C; i++ )
     for( int j = 0; j < D; j++ )
          v[i][j] = V[i][j];
          //*(v+(i*D+j)) = *(V+(i*D+j));
    return C;
}

/*=========================================================================*/

int FuzzyCMeans::WriteCentroids( char * FileName )
{
    FILE *fp;
    char  buf[1024];
    long  i, j;

    clock_t start, end;
    start = clock();

    sprintf( buf, "%s.centroids", FileName );
    printf( "\nLogging %s...", buf );
    if( (fp = fopen( buf, "w" )) == 0 )
    {
        fprintf(stderr,"%s:%d\nError opening %s for mode 'w'... [Fail]\n", __FILE__, __LINE__, buf);
        exit(1);
    }

    for( i = 0; i < C ;i++ )
    {
      for( j = 0; j < D; j++ )
           fprintf(fp, "%f\t",V[i][j]);
      fprintf(fp,"\n");
    }

    fclose(fp);

    end = clock();
    printf(" [Done in %d.%d seconds]",
    (end - start) / CLK_TCK, (end - start) % CLK_TCK);

    return 0;
}

/*=========================================================================*/

int FuzzyCMeans::WriteUmatrix( char * FileName )
{
    FILE *fp;
    char  buf[1024];
    long  i, j;

    clock_t start, end;
    start = clock();

    sprintf( buf, "%s.umatrix", FileName );
    printf( "\nLogging %s...", buf );
    if( (fp = fopen( buf, "w" )) == 0 )
    {
        fprintf(stderr,"%s:%d\nError opening %s for mode 'w'... [Fail]\n", __FILE__, __LINE__, buf);
        exit(1);
    }

    for( i = 0; i < N ;i++ )
    {
      for( j = 0; j < C; j++ )
           fprintf(fp, "%f\t",U[i][j]);
      fprintf(fp,"\n");
    }
    fclose(fp);

    end = clock();
    printf(" [Done in %d.%d seconds]",
    (end - start) / CLK_TCK, (end - start) % CLK_TCK);

    return 0;
}

/*=========================================================================*/

int FuzzyCMeans::WriteMembers( char * FileName )
{
    FILE *fp;
    char  buf[1024];
    long  i;
    int   j, max;

    clock_t start, end;
    start = clock();

    sprintf( buf, "%s.members", FileName );
    printf( "\nLogging %s...", buf );
    if( (fp = fopen( buf, "w" )) == 0 )
    {
        fprintf(stderr,"%s:%d\nError opening %s for mode 'w'... [Fail]\n", __FILE__, __LINE__, buf);
        exit(1);
    }

    for( i = 0; i < N; i++ )
    {
      for( max = j = 0; j < C; j++ ) 
         if ( U[i][j] > U[i][max] ) max = j;
            fprintf(fp,"%d\n",max);
    }
    fclose(fp);

    end = clock();
    printf(" [Done in %d.%d seconds]",
    (end - start) / CLK_TCK, (end - start) % CLK_TCK);

    return 0;
}

/*=========================================================================*/

int FuzzyCMeans::WriteClusters( char * FileName )
{
    FILE *fp;
    char  buf[1024];
    long  i, j;

    clock_t start, end;
    start = clock();

    sprintf( buf, "%s.clusters", FileName );
    printf( "\nLogging %s...", buf );
    if( (fp = fopen( buf, "w" )) == 0 )
    {
        fprintf(stderr,"%s:%d\nError opening %s for mode 'w'... [Fail]\n", __FILE__, __LINE__, buf);
        return 0;
    }

    int max;
    for( i = 0; i < N; i++ )
    {
      for( max = j = 0; j < C; j++ )
           if ( U[i][j] > U[i][max] ) max = j;
      for( j = 0; j < D; j++ )
           fprintf(fp, "%f\t", X[i][j]);
      fprintf(fp, "%d", max);
      fprintf(fp, "\n");
    }
    fclose(fp);

    end = clock();
    printf(" [Done in %d.%d seconds]",
    (end - start) / CLK_TCK, (end - start) % CLK_TCK);

    return 0;
}

/*=========================================================================*/

int FuzzyCMeans::WriteMRIClusters( int col, char * FileName )
{
    FILE *fp;
    char  buf[1024];
    long  i, j;

    clock_t start, end;
    start = clock();

    sprintf( buf, "%s.mri", FileName );
    printf( "\nLogging %s...", buf );
    if( (fp = fopen( buf, "w" )) == 0 )
    {
        fprintf(stderr,"%s:%d\nError opening %s for mode 'w'... [Fail]\n", __FILE__, __LINE__, buf);
        return 0;
    }

    int max, k = 0;
    for( i = 0; i < N; i++ )
    {
      for( max = j = 0; j < C; j++ )
           if ( U[i][j] > U[i][max] ) max = j;
      fprintf(fp, "%3d", max);
      if( k++ == (col - 1) )
      {
          fprintf(fp, "\n");
          k = 0;
      }
    }
    fclose(fp);

    end = clock();
    printf(" [Done in %d.%d seconds]",
    (end - start) / CLK_TCK, (end - start) % CLK_TCK);

    return 0;
}

/*=========================================================================*/

void FuzzyCMeans::PrintParameters( int display )
{
  if( display == 0 ) return;
  printf("\n/*=========================================================================*/");
  printf("\n/*===============================Parameters================================*/");
  printf("\n/*=========================================================================*/");
  printf("Epsilon: %f\n"    , Epsilon   );
  printf("FuzzyExp: %f\n"   , m         );
  printf("Clusters: %d\n"   , C         );
  printf("Dimension: %d\n"  , D         );
  printf("Data Length: %d\n", N         );
  printf("\n/*=========================================================================*/\n");
}

/*=========================================================================*/

void FuzzyCMeans::PrintData( int display )
{
  if( display == 0 ) return;
  printf("\n/*=========================================================================*/");
  printf("\n/*=================================Data====================================*/");
  printf("\n/*=========================================================================*/");
  for( long i = 0; i < N; i++ )
  {
     printf( "\n(%i):\t", i );
     for( int j = 0; j < D; j++ )  
        printf( "%f\t", X[i][j] );
  }
  printf("\n/*=========================================================================*/\n");
}

/***************************************************************************/
/*                      Ensemble Fuzzy C-Means Class                       */
/***************************************************************************/

/* Constructors ===========================================================*/

EnssembleFuzzyCMeans::EnssembleFuzzyCMeans()
{
        Min_Epsilon     = 0.00001;
        Max_Epsilon     = 0.00001;
        Min_m           = 2.0;
        Max_m           = 2.0;
        Min_K           = 2;
        Max_K           = 2;
        D               = 2;
        N               = 2;
        C               = 2;
        ForceC          = 0;
        FC              = 0;
        K               = 2;
        S               = 0;
        Validation      = 0;
        TC              = NULL;
        Co_Assoc        = NULL;
        X               = NULL;
        V               = NULL;
        F               = NULL;
        P               = NULL;
        SL[0] = SL[1]   = NULL;
        Consensus       = RELABELING;
        ConfusionMatrix = NULL;
        FinalClusters   = NULL;
}

EnssembleFuzzyCMeans::EnssembleFuzzyCMeans( double MinEpsilon,
                                            double MaxEpsilon,
                                            double MinFuzzyexp,
                                            double MaxFuzzyexp,
                                            int    MinCluster,
                                            int    MaxCluster,
                                            long   NrData,
                                            int    NrDimension,
                                            int    NrTrueClusters,
                                            int    FirstClusterLabel,
                                            int    NrClusterings,
                                            int    ClusterValidation,
                                            int    NrForceClusters,
                                            int    ConsensusMethod   )
{
        Min_Epsilon     = MinEpsilon;
        Max_Epsilon     = MaxEpsilon;
        Min_m           = MinFuzzyexp;
        Max_m           = MaxFuzzyexp;
        Min_K           = MinCluster;
        Max_K           = MaxCluster;
        D               = NrDimension;
        N               = NrData;
        C               = NrClusterings;
        ForceC          = NrForceClusters;
        FC              = 0;
        K               = NrTrueClusters;
        S               = FirstClusterLabel;
        Validation      = ClusterValidation;
        TC              = NULL;
        Co_Assoc        = NULL;
        X               = NULL;
        V               = NULL;
        F               = NULL;
        P               = NULL;
        SL[0]           = NULL;
        SL[1]           = NULL;
        Consensus       = ConsensusMethod;
        ConfusionMatrix = NULL;
        FinalClusters   = NULL;
}

EnssembleFuzzyCMeans::~ EnssembleFuzzyCMeans()
{
        if( V != NULL )
        {
                for( int i = 0; i < FC; i++ )
                        delete [] V[i];
                delete [] V;
                V = 0;
        }
        if( P != NULL )
        {
           for( int i = 0; i < C; i++ )
                delete [] P[i];
           delete [] P;
           P = 0;
        }
        if( FinalClusters != NULL )
        {
           for( int i = 0; i < C; i++ )
                delete [] FinalClusters[i];
           delete [] FinalClusters;
           FinalClusters = 0;
        }
        if( Co_Assoc != NULL )
        {
            delete [] Co_Assoc;
            Co_Assoc = 0;
        }
        if( F != NULL )
        {
            delete [] F;
            F = 0;
        }
        if( SL[0] != NULL )
        {
           for( long i = 0; i < N; i++ )
           {
                delete [] SL[0][i];
                delete [] SL[1][i];
           }
           delete [] SL[0];
           delete [] SL[1];
           SL[0] = SL[1] = 0;
        }
        if( Clusters != NULL )
        {
           for( int i = 0; i < 2; i++ )
                delete [] Clusters[i];
           delete [] Clusters;
           Clusters = 0;
        }
}


/* Main Methods ===========================================================*/

int EnssembleFuzzyCMeans::EFCM( char * logfile )
{
    int   i, j;
    FILE *fp = NULL;
    char  buf[512];
    sprintf( buf, "%s.log", logfile );

    clock_t start, end;
    start = clock();
    printf("\nInitialization Ensemble Fuzzy C-Means...");

    /* log all activities including time consuming */
    if( logfile != NULL )
    {
        if( (fp = fopen( buf, "a+" )) == 0 )
            fprintf(stderr,"%s:%d\nError opening %s for mode 'w'... [Fail]\n", __FILE__, __LINE__, buf);
        else
        {
            fprintf(fp, "\n>>------------------------------------------------------------------------------");
            fprintf(fp, "\nEnsemble Fuzzy C-Means on '%s' with %d samples, %d dimension, range [%d,%d] of clusters and %d clusterings", logfile, N, D, Min_K, Max_K, C);
            fprintf(fp, "\nInitialization Ensemble Fuzzy C-Means...");
        }
    }

    /* Dynamically create various data structures */
    Init();

    end = clock();
    printf(" [Done in %d.%d seconds]", (end - start) / CLK_TCK, (end - start) % CLK_TCK);

    /* log elapsed time */
    if( fp != NULL )
        fprintf(fp, " [Done in %d.%d seconds]", (end - start) / CLK_TCK, (end - start) % CLK_TCK);

    start = clock();
    printf("\nRunning Ensemble Fuzzy C-Means with %d iterations...",C);

    /* log start time */
    if( fp != NULL )
        fprintf(fp, "\nRunning Ensemble Fuzzy C-Means with %d iterations...",C);

    srand((unsigned) time(NULL));

    /* initial random cluster number for fcm between Min_K and Max_K */
    int fcm_cluster;
    if( ! ClusterGeneration )
        fcm_cluster = (rand() % (Max_K - Min_K)) + Min_K;   // randomly
    else
        fcm_cluster = Min_K;    // sequentially

    /* initial random fuzzification exponent for fcm between Min_m and Max_m */
    double fcm_fuzzyexp;
    fcm_fuzzyexp = ((int)(Max_m - Min_m) == 0)?(Min_m):((rand() % (int)(Max_m - Min_m)) + Min_m);

    /* Updates Co_Assoc iteratively */
    for( i = 0; i < C; i++ )
    {
        clock_t start2, end2;
        start2 = clock();
        printf("\n\t(%d): Running fcm with ",i);
        printf("m(%1.3f) c(%d)...",fcm_fuzzyexp,fcm_cluster);

        /* Create new fcm object by default parameters */
        fcm = new FuzzyCMeans( Min_Epsilon, fcm_fuzzyexp, fcm_cluster, D, N, 0, 0, 0, GetMaxValue() );
        fcm->SetData(X);
        fcm->FCM(0);
        if( Consensus == COASSOCIATION )
            UpdateCoAssocUmatrix(i,logfile);
        else if( Consensus == RELABELING )
            fcm->GetMembers(P[i]);       // P[i][N]
        delete fcm;
        fcm = 0;

        /* select random cluster number for fcm between Min_K and Max_K */
        if( ! ClusterGeneration )
            fcm_cluster = (rand() % (Max_K - Min_K)) + Min_K;   // randomly
        else
            if( ++fcm_cluster > Max_K ) fcm_cluster = Min_K;  // sequentially

        /* select random fuzzification exponent for fcm between Min_m and Max_m */
        fcm_fuzzyexp = ((int)(Max_m - Min_m) == 0)?(Min_m):((rand() % (int)(Max_m - Min_m)) + Min_m);

        end2 = clock();
        printf(" [Done in %d.%d seconds]", (end2 - start2) / CLK_TCK, (end2 - start2) % CLK_TCK);
    }

    end = clock();
    printf("\nRunning Ensemble Fuzzy C-Means [Done in %d.%d seconds]", (end - start) / CLK_TCK, (end - start) % CLK_TCK);

    /* log elapsed time */
    if( fp != NULL )
        fprintf(fp, " [Done in %d.%d seconds]", (end - start) / CLK_TCK, (end - start) % CLK_TCK);

    start = clock();
    if( Consensus == RELABELING )
        printf("\nExtracting Final Partition (using relabeling and voting)...");
    else if( Consensus == COASSOCIATION )
        printf("\nExtracting Final Partition (using evidence accumulation)...");

    /* log start time */
    if( fp != NULL )
        fprintf(fp, "\nExtracting Final Partition...");

    /* compute final clusters and final partition according to the consensus method */
    MakeClusters(logfile);
    
    end = clock();
    printf("\nExtracting Final Partition (%d clusters) [Done in %d.%d seconds]", FC, (end - start) / CLK_TCK, (end - start) % CLK_TCK);

    /* log elapsed time */
    if( fp != NULL )
    {
        fprintf(fp, " (%d clusters discovered)",FC);
        fprintf(fp, " [Done in %d.%d seconds]", (end - start) / CLK_TCK, (end - start) % CLK_TCK);
    }

    if( Validation )
    {
        start = clock();
        double Validity = ClusterValidation();
        printf("\nCluster Validity (Jaccard Index): %f ", Validity );
        end = clock();
        printf(" [Done in %d.%d seconds]", (end - start) / CLK_TCK, (end - start) % CLK_TCK);

        /* log cluster validity */
        if( fp != NULL )
        {
            fprintf(fp, "\nCluster Validity (Jaccard Index): %f ", Validity );
            fprintf(fp, " [Done in %d.%d seconds]", (end - start) / CLK_TCK, (end - start) % CLK_TCK);
        }
    }    

    if( fp != NULL )
        fclose( fp );

    return 0;
}

/*=========================================================================*/

/* Allocate necessary storage */
int EnssembleFuzzyCMeans::Init()
{
    if( X == NULL )
    {
        fprintf(stderr,"%s:%d\nNo storage allocated to **X... [Fail]\n", __FILE__, __LINE__);
        exit(1);
    }
        
    int i;

    /* Allocate necessary storage (Co_Assoc[N][N]) */
    if( Consensus == COASSOCIATION )
    {
        /* in this case because of we can allocate necessary storage
           from main memory, we allocate this storage from main memory */
        if( 0 )
        {
            Co_Assoc = new int [Co_AssocN];
            for( i = 0; i < Co_AssocN; i++ )
                 Co_Assoc[i] = 0;
        }
        /* in this case because of we couldn't allocate necessary storage
           from main memory, we allocate this storage from file in remainder */
        else if( 0 )
        {
            Co_Assoc = new int [N];
            for( i = 0; i < N; i++ )
                 Co_Assoc[i] = 0;
        }
        /* in this case because of we couldn't allocate necessary storage
           from main memory, we allocate this storage from file in remainder */
        else
        {
            FinalClusters = new int * [C];
            for( i = 0; i < C; i++ )
                 FinalClusters[i] = new int [N];
        }
    }

    /* Allocate necessary storage to partitions matrix(P[C][N]) */
    else if( Consensus == RELABELING )
    {
        P = new int * [C];
        for( i = 0; i < C; i++ )
             P[i] = new int [N];
    }

    return 0;
}

/*=========================================================================*/

/* Update Co_Assoc matrix */
int * EnssembleFuzzyCMeans::UpdateCoAssocUmatrix( int c, char * logfile )
{
    int * Members = new int [N];
    if( Members == NULL ) return 0;

    ////////////////////////////////////////////
    fcm->GetMembers(FinalClusters[c]); return 0;

    fcm->GetMembers(Members);

    if( 0 ) /* all computation done in main memory */
    {
        for( long i = 0; i < N - 1; i++ )
          for( int j = i + 1; j < N; j++ )
            if( Members[i] == Members[j] )
                Co_Assoc(i,j)++;
    }
    else if( 0 )   /* computation interact with file due to less memory*/
    {
        FILE *fp1,*fp2;
        char  buf1[1024];
        char  buf2[1024];

        sprintf( buf1, "%s.cooccurence", logfile );
        fp1 = fopen( buf1, "rb+" );

        sprintf( buf2, "%s.dat", "cooccurence" );
        if( (fp2 = fopen( buf2, "wb+" )) == 0 )
        {
            fprintf(stderr,"%s:%d\nError opening %s for mode 'wb+'... [Fail]\n", __FILE__, __LINE__, buf2);
            return 0;
        }

        for( long i = 0; i < N - 1; i++ )
        {
          fread( Co_Assoc, 2, 2*N, fp1 );
          for( int j = i + 1; j < N; j++ )
            if( Members[i] == Members[j] )
                Co_Assoc[j]++;
          fwrite( Co_Assoc, 2, 2*N, fp2 );
        }
        fclose(fp1);
        fclose(fp2);
        remove(buf1);
        rename(buf2,buf1);
    }
    else
    {
        for( long i = 0; i < N - 1; i++ )
             FinalClusters[c][i] = Members[i];
    }

    delete [] Members;
    Members = 0;
    return Co_Assoc;
}

/*=========================================================================*/

int EnssembleFuzzyCMeans::MakeClusters( char * logfile )
{
    F  = new int [N];  /* final cluster label for every data */
    FC = 0;            /* number of final clusters */

    long i, j;

    if( Consensus == RELABELING )
    {
        for( i = 0; i < N; i++ )
        {
          int maxcount = MININT;
          for( j = 0; j < C; j++ )
          {
             int count = 0;
             for( int k = 0; k < C; k++ )
               if( P[k][i] == P[j][i] )
                   count++;
             if( count > maxcount )
             {
               F[i] = P[j][i];
               maxcount = count;
             }
          }
          /* final cluster number (in fact: FC is the maximum number in **P) */
          if( F[i] > FC )
              FC = F[i];
        }
        FC++;
    }
    else if( Consensus == COASSOCIATION )
    {
        /* Allocate necessary storage to SL[2][N][C] (Single-Linkage-Tree) */
        SL[0] = new long * [N];
        SL[1] = new long * [N];
        for( i = 0; i < N; i++ )
        {
             SL[0][i] = new long [C+1];
             SL[1][i] = new long [C+1];
        }

        for( i = 0; i < C + 1; i++ )
          for( j = 0; j < N; j++ )
          {
             SL[0][j][i] = UNKNOWNINT;
             SL[1][j][i] = UNKNOWNINT;
          }

        /* Allocate necessary storage to **Clusters */
        Clusters = new long * [2];
        for( i = 0; i < 2; i++ ) 
             Clusters[i] = new long [C+1];

        /* long Clusters[2][C+1]; 
           store number of clusters and level of each partition. 
           first dimension contains two values involve 
           cluster-number ([0]) and level of partition ([1]).
           second dimension is the maximum partition number ([C+1]),
           actually is the number of ensembles */

        Clustersidx      = 1;  /* index of Clusters[2][C+1] (used in second dimension) */
        Clusters[0][0]   = N;  /* at first, number of clusters is N (length of data) */
        Clusters[1][0]   = 0;  /* at first, level of first partition is 0 */
        int Level        = 0;  /* current level of SL (single-linkage-tree) */

        /* create single-linkage tree on co-association matrix */
        printf("\n   Constructing Single-Linkage-Tree with %d levels...",C + 1);
        while( Level < C + 1 )
        {
          clock_t start, end;
          start = clock();
          printf("\n\t(%d):",Level);

          if( 0 ) /* all computation done in main memory */
          {
           for( i = 0; i < N - 1; i++ )
            for( j = i + 1; j < N; j++ )
               if( Co_Assoc(i,j) == (C - Level) )
               {
                   int t1 = i, t2 = j;
                   /* go back to first element in which cluster that contain i */
                   while( SL[0][t1][Level] != UNKNOWNINT )
                          t1 = SL[0][t1][Level];
                   /* go back to first element in which cluster that contain j */
                   while( SL[0][t2][Level] != UNKNOWNINT )
                          t2 = SL[0][t2][Level];

                   /* merge two clusters that contain i and j */
                   while( t1 != UNKNOWNINT || t2 != UNKNOWNINT )
                   {
                          /* data i, j are already in the same cluster */
                          if( t1 == t2 )
                              break;
                          else if( t1 < t2 )
                          {
                              while( SL[1][t1][Level] != UNKNOWNINT &&
                                     SL[1][t1][Level] < t2 )
                                     t1 = SL[1][t1][Level];
                              int t = SL[1][t1][Level];
                              SL[1][t1][Level] = t2;
                              SL[0][t2][Level] = t1;
                              if( t == UNKNOWNINT ) 
                                  break;
                              t1 = t;
                          }
                          else if( t1 > t2 )
                          {
                              while( SL[1][t2][Level] != UNKNOWNINT &&
                                     SL[1][t2][Level] < t1 )
                                     t2 = SL[1][t2][Level];
                              int t = SL[1][t2][Level];
                              SL[1][t2][Level] = t1;
                              SL[0][t1][Level] = t2;
                              if( t == UNKNOWNINT ) 
                                  break;
                              t2 = t;
                          }
                   } // while
               } // if
          }
          else if( 0 )   /* computation interact with file due to less memory*/
          {
               FILE *fp;
               char  buf[1024];

               sprintf( buf, "%s.cooccurence", logfile );
               if( (fp = fopen( buf, "rb+" )) == 0 )
               {
                   fprintf(stderr,"%s:%d\nError opening %s for mode 'wb+'... [Fail]\n", __FILE__, __LINE__, buf);
                   //return 0;
               }

               for( i = 0; i < N - 1; i++ )
               {
                 fread( Co_Assoc, 2, 2*N, fp );
                 for( j = i + 1; j < N; j++ )
                 if( Co_Assoc[j] == (C - Level) )
                 {
                   int t1 = i, t2 = j;
                   /* go back to first element in which cluster that contain i */
                   while( SL[0][t1][Level] != UNKNOWNINT )
                          t1 = SL[0][t1][Level];
                   /* go back to first element in which cluster that contain j */
                   while( SL[0][t2][Level] != UNKNOWNINT )
                          t2 = SL[0][t2][Level];

                   /* merge two clusters that contain i and j */
                   while( t1 != UNKNOWNINT || t2 != UNKNOWNINT )
                   {
                          /* data i, j are already in the same cluster */
                          if( t1 == t2 )
                              break;
                          else if( t1 < t2 )
                          {
                              while( SL[1][t1][Level] != UNKNOWNINT &&
                                     SL[1][t1][Level] < t2 )
                                     t1 = SL[1][t1][Level];
                              int t = SL[1][t1][Level];
                              SL[1][t1][Level] = t2;
                              SL[0][t2][Level] = t1;
                              if( t == UNKNOWNINT ) 
                                  break;
                              t1 = t;
                          }
                          else if( t1 > t2 )
                          {
                              while( SL[1][t2][Level] != UNKNOWNINT &&
                                     SL[1][t2][Level] < t1 )
                                     t2 = SL[1][t2][Level];
                              int t = SL[1][t2][Level];
                              SL[1][t2][Level] = t1;
                              SL[0][t1][Level] = t2;
                              if( t == UNKNOWNINT ) 
                                  break;
                              t2 = t;
                          }
                   } // while
                 } // if
               }
               fclose(fp);
          }
          else            // FinalClusters[C][N]
          {
           Co_Assoc = new int [N];
           for( i = 0; i < N - 1; i++ )
           {
                for( j = i + 1; j < N; j++ ) Co_Assoc[j] = 0;
                for( int c = 0; c < C; c++ )
                for( j = i + 1; j < N; j++ )
                if( FinalClusters[c][i] == FinalClusters[c][j] )
                Co_Assoc[j]++;
            for( j = i + 1; j < N; j++ )
               if( Co_Assoc[j] == (C - Level) )
               {
                   int t1 = i, t2 = j;
                   /* go back to first element in which cluster that contain i */
                   while( SL[0][t1][Level] != UNKNOWNINT )
                          t1 = SL[0][t1][Level];
                   /* go back to first element in which cluster that contain j */
                   while( SL[0][t2][Level] != UNKNOWNINT )
                          t2 = SL[0][t2][Level];

                   /* merge two clusters that contain i and j */
                   while( t1 != UNKNOWNINT || t2 != UNKNOWNINT )
                   {
                          /* data i, j are already in the same cluster */
                          if( t1 == t2 )
                              break;
                          else if( t1 < t2 )
                          {
                              while( SL[1][t1][Level] != UNKNOWNINT &&
                                     SL[1][t1][Level] < t2 )
                                     t1 = SL[1][t1][Level];
                              int t = SL[1][t1][Level];
                              SL[1][t1][Level] = t2;
                              SL[0][t2][Level] = t1;
                              if( t == UNKNOWNINT ) 
                                  break;
                              t1 = t;
                          }
                          else if( t1 > t2 )
                          {
                              while( SL[1][t2][Level] != UNKNOWNINT &&
                                     SL[1][t2][Level] < t1 )
                                     t2 = SL[1][t2][Level];
                              int t = SL[1][t2][Level];
                              SL[1][t2][Level] = t1;
                              SL[0][t1][Level] = t2;
                              if( t == UNKNOWNINT ) 
                                  break;
                              t2 = t;
                          }
                   } // while
               } // if
           }
           delete [] Co_Assoc;
          }

          Level++;
          Clusters[0][Clustersidx] = 0;
          if( Level < C + 1 )
          {
            for( i = 0; i < N; i++ )
            {
                 /* to initiate current partition, transfer previous 
                    partition (Level-1) to current partition (Level) */
                 SL[1][i][Level] = SL[1][i][Level-1];
                 SL[0][i][Level] = SL[0][i][Level-1];

                 /* find cluster number for each level of SL-tree */
                 if( SL[1][i][Level] == UNKNOWNINT )
                     Clusters[0][Clustersidx]++;
            }
          }
          else
          {
              Clusters[0][Clustersidx] = 1;
              Clusters[1][Clustersidx] = Level;
          }

          if( Clusters[0][Clustersidx] < Clusters[0][Clustersidx-1] )
              Clusters[1][Clustersidx++] = Level - 1;

          end = clock();
          printf(" [Done in %d.%d seconds]", (end - start) / CLK_TCK, (end - start) % CLK_TCK);
        } // while( Level < C + 1 )

        /* find maximum lifetime in single-linkage-tree */
        int lifetime = 0;
        int lifetimelevel = 0;
        for( i = 0; i < Clustersidx - 1; i++ )
           if( Clusters[1][i+1] - Clusters[1][i] >= lifetime )
           {
               lifetime = Clusters[1][i+1] - Clusters[1][i];
               lifetimelevel = Clusters[1][i];
               FC = Clusters[0][i];
           }

        /* force to produce specified clusters */
        if( ForceC > 0 )
        for( i = 0; i < Clustersidx; i++ )
           if( Clusters[0][i] == ForceC )
           {
               lifetime = Clusters[1][i+1] - Clusters[1][i];
               lifetimelevel = Clusters[1][i];
               FC = Clusters[0][i];
           }

        /* put final cluster label for each data (final partition) */
        for( i = 0; i < N; i++ )
             F[i] = UNKNOWNINT;
        int clusterlabel = 0;
        for( i = 0; i < N; i++ )
             if( F[i] == UNKNOWNINT )
             {
                 int t = i;
                 do{
                        F[t] = clusterlabel;
                        t = SL[1][t][lifetimelevel];
                 }while( t != UNKNOWNINT );
                 clusterlabel++;
             }
    }

    /* make cluster prototypes */
    /* Allocate necessary storage to **V (centroids) */
    if( V == NULL )
    {
        V = new double * [FC];
        for( i = 0; i < FC; i++ )
             V[i] = new double [D];
    }

    for( i = 0; i < FC; i++ )
     for( j = 0; j < D; j++ )
          V[i][j] = 0;

    for( int k = 0; k < FC; k++ )
    { 
      for( i = 0; i < D; i++ )
      {
           int Denominator = 0;
           for( j = 0; j < N; j++ )
           if( F[j] == k )
           {
               V[k][i] += X[j][i];
               Denominator++;
           }
           (Denominator > 0) ? (V[k][i] /= Denominator) : (V[k][i] = MAX);
      }
    }

    return 0;
}

/*=========================================================================*/

double EnssembleFuzzyCMeans::ClusterValidation( )
{
    int i, j;

    /* Allocate necessary storage to **Clusters */
    ConfusionMatrix = new long * [K];
    for( i = 0; i < K; i++ ) 
         ConfusionMatrix[i] = new long [FC];

    for( i = 0; i < K; i++ )
      for( j = 0; j < FC; j++ ) 
         ConfusionMatrix[i][j] = 0;

    for( i = 0; i < K; i++ )
      for( j = 0; j < FC; j++ )
      {
           for( long x = 0; x < N; x++ )
             if( TC[x] == i + S )
             {
               for( long y = 0; y < N; y++ )
                 if( F[y] == j )
                   if( x == y )
                     ConfusionMatrix[i][j]++;
             }
      }

    long double N11, N10, N01, N00;
    /* The category N11 contains the pairs of points that are in the same 
       cluster both in C and in C'. The category N10 contains the pairs of 
       points that are in the same cluster in C but not in C'. The 
       definitions of N01 and N00 are similar. */

    N11 = N10 = N01 = N00 = 0;

    // ----------N11
    for( i = 0; i < K; i++ )
      for( j = 0; j < FC; j++ )
        N11 += ConfusionMatrix[i][j] * ConfusionMatrix[i][j];
    N11 -= N;
    N11 /= 2;

    // ----------N10
    long n = 0;
    for( j = 0; j < FC; j++ )
    {
      n = 0;
      for( i = 0; i < K; i++ )
         n += ConfusionMatrix[i][j];
      N10 += n * n;
    }

    n = 0;
    for( i = 0; i < K; i++ )
      for( j = 0; j < FC; j++ )
         n += ConfusionMatrix[i][j] * ConfusionMatrix[i][j];

    N10 -= n;
    N10 /= 2;

    // ----------N01
    n = 0;
    for( i = 0; i < K; i++ )
    {
      n = 0;
      for( j = 0; j < FC; j++ )
         n += ConfusionMatrix[i][j];
      N01 += n * n;
    }

    n = 0;
    for( i = 0; i < K; i++ )
      for( j = 0; j < FC; j++ )
         n += ConfusionMatrix[i][j] * ConfusionMatrix[i][j];

    N01 -= n;
    N01 /= 2;

    // ----------N00
    long double PP;
    double Jaccard;
    /* PP = N(N - 1)/2 is the total number of point pairs */

    PP = N * ( N - 1 ) / 2;
    N00 = PP - N11 - N10 - N01;

    Jaccard = N11 / ( N11 + N01 + N10 );

    for( i = 0; i < K; i++ )
         delete [] ConfusionMatrix[i];
    delete [] ConfusionMatrix;
    ConfusionMatrix = 0;

    return Jaccard;
}

/*=========================================================================*/

int EnssembleFuzzyCMeans::GetNrFinalCluster()
{
    return FC;
}

/*=========================================================================*/

int EnssembleFuzzyCMeans::GetFinalClusters( int * f )
{
    for( long i = 0; i < N ;i++ )
         *(f+i) = *(F+i);
    return FC;
}

/*=========================================================================*/

int EnssembleFuzzyCMeans::GetClusters( double **v )
{
    for( int i = 0; i < FC; i++ )
     for( int j = 0; j < D; j++ )
          v[i][j] = V[i][j];
          //*(v+(i*D+j)) = *(V+(i*D+j));
    return FC;
}

/* Writing Out Methods ====================================================*/

int EnssembleFuzzyCMeans::WriteCoAssocMatrix( char * FileName )
{
    clock_t start, end;
    start = clock();

    FILE *fp;
    char  buf[1024];
    long  i, j;

    sprintf( buf, "%s.coassoc", FileName );
    printf( "\nLogging %s...", buf );
    if( (fp = fopen( buf, "w" )) == 0 )
    {
        fprintf(stderr,"%s:%d\nError opening %s for mode 'w'... [Fail]\n", __FILE__, __LINE__, buf);
        return 0;
    }

    for( i = 0; i < N ;i++ )
    {
      for( j = 0; j < N; j++ )
      {
         if( i == j )
             fprintf(fp, "-\t");
         else if( i > j )
             fprintf(fp, "%1.2f\t", (float)Co_Assoc(j,i)/C);
         else
             fprintf(fp, "%1.2f\t", (float)Co_Assoc(i,j)/C);
      }
      fprintf(fp,"\n");
    }
    fclose(fp);

    end = clock();
    printf(" [Done in %d.%d seconds]",
    (end - start) / CLK_TCK, (end - start) % CLK_TCK);
  
    return 0;
}

/*=========================================================================*/

int EnssembleFuzzyCMeans::WriteSLTree( char * FileName )
{
    FILE *fp;
    char  buf[1024];
    long  i, j;

    clock_t start, end;
    start = clock();
      
    sprintf( buf, "%s.sltree", FileName );
    printf( "\nLogging %s...", buf );
    if( (fp = fopen( buf, "w" )) == 0 )
    {
        fprintf(stderr,"%s:%d\nError opening %s for mode 'w'... [Fail]\n", __FILE__, __LINE__, buf);
        return 0;
    }

    /* write sltree from top to bottom, in which each level of sltree writed in
       each row of file */
    for( i = 0; i < C + 1; i++ )
    {
      for( j = 0; j < N; j++ )
         if( SL[1][j][i] != UNKNOWNINT )
              fprintf(fp, "%5d", SL[1][j][i]+1);
         else
              fprintf(fp, "%5d", -1);
        fprintf(fp, "\n");
    }

    fclose(fp);
    end = clock();
    printf(" [Done in %d.%d seconds]",
    (end - start) / CLK_TCK, (end - start) % CLK_TCK);

    start = clock();

    sprintf( buf, "%s.sltreelifetimes", FileName );
    printf( "\nLogging %s...", buf );
    if( (fp = fopen( buf, "w" )) == 0 )
    {
        fprintf(stderr,"%s:%d\nError opening %s for mode 'w'... [Fail]\n", __FILE__, __LINE__, buf);
        return 0;
    }

    fprintf(fp, "Level of SL-Tree:\t");
    for( j = 0; j < Clustersidx; j++ )
         fprintf(fp, "%4d", Clusters[1][j]);
    fprintf(fp, "\nClusters Number:\t");
    for( j = 0; j < Clustersidx; j++ )
         fprintf(fp, "%4d", Clusters[0][j]);
    fprintf(fp, "\n");

    fclose(fp);
    end = clock();
    printf(" [Done in %d.%d seconds]",
    (end - start) / CLK_TCK, (end - start) % CLK_TCK);

    return 0;
}

/*=========================================================================*/

int EnssembleFuzzyCMeans::WriteMembers( char * FileName )
{
    FILE *fp;
    char  buf[1024];
    long  i;

    clock_t start, end;
    start = clock();
      
    sprintf( buf, "%s.members", FileName );
    printf( "\nLogging %s...", buf );
    if( (fp = fopen( buf, "w" )) == 0 )
    {
        fprintf(stderr,"%s:%d\nError opening %s for mode 'w'... [Fail]\n", __FILE__, __LINE__, buf);
        return 0;
    }

    for( i = 0; i < N; i++ )
         fprintf(fp, "%d\n", F[i]);

    fclose(fp);
    end = clock();
    printf(" [Done in %d.%d seconds]",
    (end - start) / CLK_TCK, (end - start) % CLK_TCK);

    return 0;
}

/*=========================================================================*/

int EnssembleFuzzyCMeans::WriteCentroids( char * FileName )
{
    FILE *fp;
    char  buf[1024];
    int   i, k;
    long  j;

    clock_t start, end;
    start = clock();
      
    sprintf( buf, "%s.centroids", FileName );
    printf( "\nLogging %s...", buf );
    if( (fp = fopen( buf, "w" )) == 0 )
    {
        fprintf(stderr,"%s:%d\nError opening %s for mode 'w'... [Fail]\n", __FILE__, __LINE__, buf);
        return 0;
    }

    for( i = 0; i < FC; i++ )
    {
      for( j = 0; j < D; j++ )
      {
           if( V[i][j] != MAX ) fprintf(fp, "%f\t", V[i][j]);
      }
      fprintf(fp, "\n");
    }

    fclose(fp);
    end = clock();
    printf(" [Done in %d.%d seconds]",
    (end - start) / CLK_TCK, (end - start) % CLK_TCK);

    return 0;
}

/*=========================================================================*/

int EnssembleFuzzyCMeans::WriteClusters( char * FileName, int integration )
{
    FILE *fp;
    char  buf[1024];
    long  i;
    int   j, k;

    if( integration == JOINED )
    {
           clock_t start, end;
           start = clock();

           sprintf( buf, "%s.clusters", FileName );
           printf( "\nLogging %s...", buf );
           if( (fp = fopen( buf, "w" )) == 0 )
           {
               fprintf(stderr,"%s:%d\nError opening %s for mode 'w'... [Fail]\n", __FILE__, __LINE__, buf);
               return 0;
           }

           for( i = 0; i < N; i++ )
           {
                for( j = 0; j < D; j++ )
                     fprintf(fp, "%f\t", X[i][j]);
                fprintf(fp, "%d", F[i]);
                fprintf(fp, "\n");
           }
           fclose(fp);

           end = clock();
           printf(" [Done in %d.%d seconds]",
           (end - start) / CLK_TCK, (end - start) % CLK_TCK);
    }
    else // if( integration == DISJOINED )
    {
           char *cluster[] = {"1","2","3","4","5","6","7","8","9","10",
                              "11","12","13","14","15","16","17","18","19","20",
                              "21","22","23","24","25","26","27","28","29","30"};

           for( k = 0; k < FC; k++ )
           {
                // check does cluster k contain any data
                int b = 1;
                for( i = 0; i < N; i++ )
                     if( F[i] == k )
                     {
                         b = 0;
                         break;
                     }
                if( b ) continue;

                clock_t start, end;
                start = clock();

                sprintf( buf, "%s.cluster", FileName );
                //char num[10];
                //itoa(k+1,num,10);
                //strcat(buf,num);
                strcat(buf,cluster[k]);
                printf( "\nLogging %s...", buf );
                if( (fp = fopen( buf, "w" )) == 0 )
                {
                    fprintf(stderr,"%s:%d\nError opening %s for mode 'w'... [Fail]\n", __FILE__, __LINE__, buf);
                    return 0;
                }

                for( i = 0; i < N; i++ )
                     if( F[i] == k )
                     {
                         for( j = 0; j < D; j++ )
                              fprintf(fp, "%f\t", X[i][j]);
                         fprintf(fp, "\n");
                     }
                fclose(fp);

                end = clock();
                printf(" [Done in %d.%d seconds]",
                (end - start) / CLK_TCK, (end - start) % CLK_TCK);
           }
    }

    return 0;
}

/*=========================================================================*/

int EnssembleFuzzyCMeans::WriteMRIClusters( int col, char * FileName )
{
    FILE *fp;
    char  buf[1024];
    long  i;

    clock_t start, end;
    start = clock();

    sprintf( buf, "%s.mri", FileName );
    printf( "\nLogging %s...", buf );
    if( (fp = fopen( buf, "w" )) == 0 )
    {
        fprintf(stderr,"%s:%d\nError opening %s for mode 'w'... [Fail]\n", __FILE__, __LINE__, buf);
        return 0;
    }

    int k = 0;
    for( i = 0; i < N; i++ )
    {
      fprintf(fp, "%3d", F[i]);
      if( k++ == (col - 1))
      {
          fprintf(fp, "\n");
          k = 0;
      }
    }
    fclose(fp);

    end = clock();
    printf(" [Done in %d.%d seconds]",
    (end - start) / CLK_TCK, (end - start) % CLK_TCK);

    return 0;
}

/*=========================================================================*/

int EnssembleFuzzyCMeans::WritePartitions( char * FileName )
{
    if( P == NULL ) return 1;
    FILE *fp;
    char  buf[1024];
    long  i;
    int   j, k;

    clock_t start, end;
    start = clock();

    sprintf( buf, "%s.partitions", FileName );
    printf( "\nLogging %s...", buf );
    if( (fp = fopen( buf, "w" )) == 0 )
    {
        fprintf(stderr,"%s:%d\nError opening %s for mode 'w'... [Fail]\n", __FILE__, __LINE__, buf);
        return 0;
    }

    for( i = 0; i < N; i++ )
    {
       for( j = 0; j < C; j++ )
                fprintf(fp, "%3d", P[j][i]);
       fprintf(fp, "\n");
    }
    fclose(fp);

    end = clock();
    printf(" [Done in %d.%d seconds]",
    (end - start) / CLK_TCK, (end - start) % CLK_TCK);

    return 0;
}

/* Loading Data ===========================================================*/

void EnssembleFuzzyCMeans::SetData( double **x )
{
        X = x;
}

/*=========================================================================*/

void EnssembleFuzzyCMeans::SetTrueClusters( int * tc )
{
        TC = tc;
}

/*=========================================================================*/

void EnssembleFuzzyCMeans::SetMaxValue( double * maxvalue )
{
        MaxValue = maxvalue;
}

/*=========================================================================*/

double * EnssembleFuzzyCMeans::GetMaxValue()
{
        return MaxValue;
}

/*=========================================================================*/

void EnssembleFuzzyCMeans::PrintParameters( int display )
{
  if( display == 0 ) return;
  printf("\n/*=================Ensemble Clustering Parameters===================*/\n");
  printf("[MinimumEpsilon]:\t%f\n" , Min_Epsilon );
  printf("[MaximumEpsilon]:\t%f\n" , Max_Epsilon );
  printf("[MinimumFuzzyExp]:\t%f\n", Min_m       );
  printf("[MaximumFuzzyExp]:\t%f\n", Max_m       );
  printf("[MinimumClusters]:\t%d\n", Min_K       );
  printf("[MaximumClusters]:\t%d\n", Max_K       );
  printf("[DataLength]:\t\t%d\n"   , N           );
  printf("[Dimensions]:\t\t%d\n"   , D           );
  printf("[Clusterings]:\t\t%d\n"  , C           );
  printf("[ConsensusMethod]:\t%d\n", Consensus   );
  printf("/*====================================================================*/\n");
}

/*=========================================================================*/

/***************************************************************************/
/*                  Stream Ensemble Fuzzy C-Means Class                    */
/***************************************************************************/

/* Constructors ===========================================================*/

StreamEnssembleFuzzyCMeans::StreamEnssembleFuzzyCMeans()
{
        Min_Epsilon     = 0.00001;
        Max_Epsilon     = 0.00001;
        Min_m           = 2.0;
        Max_m           = 2.0;
        Min_K           = 2;
        Max_K           = 2;
        D               = 2;
        N               = 2;
        B               = 1000;
        C               = 2;
        ForceC          = 0;
        FC              = 0;
        K               = 2;
        S               = 0;
        Validation      = 0;
        TC              = NULL;
        X               = NULL;
        V               = NULL;
        F               = NULL;
        ConfusionMatrix = NULL;
        Consensus       = RELABELING;
        BaseMethod      = new char [10];
        strcpy( BaseMethod, _EFCM );
        BIdx            = 0;
        NrBlocks        = 0;
        BFC             = NULL;
        BF              = NULL;
        BV              = NULL;
        BL              = NULL;
        BTC             = NULL;
}

StreamEnssembleFuzzyCMeans::StreamEnssembleFuzzyCMeans( double MinEpsilon,
                                                        double MaxEpsilon,
                                                        double MinFuzzyexp,
                                                        double MaxFuzzyexp,
                                                        int    MinCluster,
                                                        int    MaxCluster,
                                                        long   NrData,
                                                        int    NrDimension,
                                                        long   BlockSize,
                                                        int    NrTrueClusters,
                                                        int    FirstClusterLabel,
                                                        int    NrClusterings,
                                                        int    ClusterValidation,
                                                        int    NrForceClusters,
                                                        int    ConsensusMethod,
                                                        char * BaseClusterMethod )
{
        Min_Epsilon     = MinEpsilon;
        Max_Epsilon     = MaxEpsilon;
        Min_m           = MinFuzzyexp;
        Max_m           = MaxFuzzyexp;
        Min_K           = MinCluster;
        Max_K           = MaxCluster;
        D               = NrDimension;
        N               = NrData;
        B               = BlockSize;
        C               = NrClusterings;
        ForceC          = NrForceClusters;
        FC              = 0;
        K               = NrTrueClusters;
        S               = FirstClusterLabel;
        Validation      = ClusterValidation;
        TC              = NULL;
        X               = NULL;
        V               = NULL;
        F               = NULL;
        ConfusionMatrix = NULL;
        Consensus       = ConsensusMethod;
        BaseMethod      = new char [10];
        strcpy( BaseMethod, BaseClusterMethod );
        BIdx            = 0;
        NrBlocks        = (N % B == 0)?(N / B):((N / B) + 1);
        BFC             = new int       [NrBlocks];
        BF              = new int    *  [NrBlocks];
        BV              = new double ** [NrBlocks];
        BL              = new long      [NrBlocks];
        BTC             = new int    *  [NrBlocks];
        for(int i=0;i<NrBlocks;i++) BTC[i] = new int [B];
}

StreamEnssembleFuzzyCMeans::~ StreamEnssembleFuzzyCMeans()
{
    /*
        if( V != NULL )
        {
            for( int i = 0; i < FC; i++ )
                 delete [] V[i];
            delete [] V;
            V = 0;
        }
        if( F != NULL )
        {
            delete [] F;
            F = 0;
        }
        if( X != NULL )
        {
            for( int i = 0; i < D; i++ )
                 delete [] X[i];
            delete [] X;
            X = NULL;
        }
        if( TC != NULL )
        {
            delete [] TC;
            TC = 0;
        }
    */
}

/* Main Methods ===========================================================*/

int StreamEnssembleFuzzyCMeans::Init( char * logfile )
{
    FILE *fp = NULL;
    char  buf[512];
    sprintf( buf, "%s.log", logfile );

    clock_t start, end;
    start = clock();
    printf("\nInitialization Stream Ensemble Fuzzy C-Means...");

    /* log all activities including time consuming */
    if( logfile != NULL )
    {
        if( (fp = fopen( buf, "a+" )) == 0 )
            fprintf(stderr,"%s:%d\nError opening %s for mode 'w'... [Fail]\n", __FILE__, __LINE__, buf);
        else
        {
            fprintf(fp, "\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
            fprintf(fp, "\nStream Ensemble Fuzzy C-Means on '%s' with %d samples, %d dimension, block size %d, range [%d,%d] of clusters and %d clusterings", logfile, N, D, B, Min_K, Max_K, C);
            fprintf(fp, "\nInitialization Stream Ensemble Fuzzy C-Means...");
        }
    }

    end = clock();
    printf(" [Done in %d.%d seconds]", (end - start) / CLK_TCK, (end - start) % CLK_TCK);

    /* log elapsed time */
    if( fp != NULL )
        fprintf(fp, " [Done in %d.%d seconds]", (end - start) / CLK_TCK, (end - start) % CLK_TCK);

    fclose(fp);
    return 0;
}

/*=========================================================================*/

int StreamEnssembleFuzzyCMeans::SEFCM( long Block, char * logfile )
{
    long  i, j;
    FILE *fp = NULL;
    char  buf[512];
    sprintf( buf, "%s.log", logfile );

    clock_t start, end;
    start = clock();

    printf("\nRunning Stream Ensemble Fuzzy C-Means on Block %d with %d samples...",BIdx, Block);

    /* log all activities including time consuming */
    if( logfile != NULL )
    {
        if( (fp = fopen( buf, "a+" )) == 0 )
            fprintf(stderr,"%s:%d\nError opening %s for mode 'w'... [Fail]\n", __FILE__, __LINE__, buf);
    }

    if( !strcmp(BaseMethod, "efcm") )
    {
        /* Create new efcm object by default parameters */
        EnssembleFuzzyCMeans * efcm;
        efcm = new EnssembleFuzzyCMeans( Min_Epsilon, Max_Epsilon, Min_m, Max_m, Min_K, Max_K, 
                                         Block, D, K, S, C, Validation, ForceC, Consensus );

        GetTrueClusters      ( BTC[BIdx]     );

        /* Set data */
        efcm->SetData        ( GetData()     );
        efcm->SetTrueClusters( BTC[BIdx]     );
        efcm->SetMaxValue    ( GetMaxValue() );

        /* Run ensemble fuzzy c-means clustering */
        efcm->EFCM           ( logfile       );

        /* store result of current block clustering */
        BL [BIdx] = Block;
        BFC[BIdx] = efcm->GetNrFinalCluster();
        BF [BIdx] = new int [Block];
        BV [BIdx] = new double * [BFC[BIdx]];
        for( j = 0; j < BFC[BIdx]; j++ )
             BV[BIdx][j] = new double [D];

        efcm->GetFinalClusters( BF [BIdx] );
        efcm->GetClusters     ( BV [BIdx] );

        delete efcm;
    }
    else // if( !strcmp(BaseMethod, "fcm") )
    {

        /* initial random cluster number for fcm between Min_K and Max_K */
        srand((unsigned) time(NULL));
        static int fcm_cluster = Min_K;    // sequentially

        /* select random cluster number for fcm between Min_K and Max_K */
        if( ++fcm_cluster > Max_K ) fcm_cluster = Min_K;  // sequentially

        /* Create new fcm object by default parameters */
        FuzzyCMeans * fcm;
        fcm = new FuzzyCMeans( Min_Epsilon, Min_m, ((rand() % (Max_K - Min_K)) + Min_K), D,
                               Block, K, Validation, S );

        GetTrueClusters      ( BTC[BIdx]     );

        /* Set data */
        fcm->SetData         ( GetData()     );
        fcm->SetTrueClusters ( BTC[BIdx]     );
        fcm->SetMaxValue     ( GetMaxValue() );

        /* Run fuzzy c-means clustering */
        clock_t start2, end2;
        start2 = clock();
        printf("\nRunning Fuzzy C-Means...");

        /* log all activities including time consuming */
        if( fp != NULL )
        {
            fprintf(fp, "\n>>------------------------------------------------------------------------------");
            fprintf(fp, "\nRunning Fuzzy C-Means on '%s' with %d samples, %d dimension and %d clusters",logfile, Block, D, Min_K);
        }

        /* Run fuzzy c-means clustering */
        fcm->FCM(1);

        end2 = clock();
        printf(" [Done in %d.%d seconds]", (end2 - start2) / CLK_TCK, (end2 - start2) % CLK_TCK);

        /* log elapsed time */
        if( fp != NULL )
            fprintf(fp, " [Done in %d.%d seconds]", (end2 - start2) / CLK_TCK, (end2 - start2) % CLK_TCK);

        if( Validation )
        {
            start2 = clock();
            double Validity = fcm->ClusterValidation();
            printf("\nCluster Validity (Jaccard Index): %f ", Validity );
            end2   = clock();
            printf(" [Done in %d.%d seconds]", (end2 - start2) / CLK_TCK, (end2 - start2) % CLK_TCK);

            /* log cluster validity */
            if( fp != NULL )
            {
                fprintf(fp, "\nCluster Validity (Jaccard Index): %f ", Validity );
                fprintf(fp, " [Done in %d.%d seconds]", (end2 - start2) / CLK_TCK, (end2 - start2) % CLK_TCK);
            }
        }

        /* store result of current block clustering */
        BL [BIdx] = Block;
        BFC[BIdx] = fcm->GetNrFinalCluster();
        BF [BIdx] = new int [Block];
        BV [BIdx] = new double * [BFC[BIdx]];
        for( j = 0; j < BFC[BIdx]; j++ )
             BV[BIdx][j] = new double [D];

        fcm->GetFinalClusters( BF [BIdx] );
        fcm->GetClusters     ( BV [BIdx] );

        delete fcm;
    }

    BIdx++;

    end = clock();
    printf("\nRunning SEFCM on Block %d with %d samples...",BIdx, Block);
    printf(" [Done in %d.%d seconds]", (end - start) / CLK_TCK, (end - start) % CLK_TCK);

    /* log elapsed time */
    if( fp != NULL )
    {
        fprintf(fp,"\nRunning Stream Ensemble Fuzzy C-Means on Block %d with %d samples...",BIdx, Block);
        fprintf(fp, " [Done in %d.%d seconds]", (end - start) / CLK_TCK, (end - start) % CLK_TCK);
    }

    fclose( fp );

    return 0;
}

/*=========================================================================*/

int StreamEnssembleFuzzyCMeans::ReCluster( char * logfile )
{
    long  i, j, k, x, y;
    FILE *fp = NULL;
    char  buf[512];
    sprintf( buf, "%s.log", logfile );

    clock_t start, end;
    start = clock();

    /* sum of all clusters at every iteration (as input data for reclustering) */
    int TotalClusters = 0;
    for( i = 0; i < BIdx; i++ )
         TotalClusters += BFC[i];

    printf("\nReclustering result of %d blocks(iterations) with %d samples(clusters)...",BIdx, TotalClusters);

    /* log all activities including time consuming */
    if( logfile != NULL )
    {
        if( (fp = fopen( buf, "a+" )) == 0 )
            fprintf(stderr,"%s:%d\nError opening %s for mode 'w'... [Fail]\n", __FILE__, __LINE__, buf);
        else
            fprintf(fp,"\nReclustering result of %d blocks(iterations) with %d samples(clusters)...",BIdx, TotalClusters);
    }

    /* allocating necessary storage to X for reclustering */
    X = new double * [TotalClusters];
    for( i = 0; i < TotalClusters; i++ )
         X[i] = new double [D];

    /* casting cluster centers (***BV) into **X as input data for reclustering */
    for( x = 0, k =0; x < BIdx; x++ )
      for( i = 0; i < BFC[x]; i++, k++ )
        for( j = 0; j < D; j++ )
             X[k][j] = BV[x][i][j];

    /* Create new efcm object by default parameters */
    EnssembleFuzzyCMeans * efcm;
    efcm = new EnssembleFuzzyCMeans( Min_Epsilon, Max_Epsilon, Min_m, Max_m, Min_K, Max_K, 
                                     TotalClusters, D, K, S, C, 0, ForceC, Consensus );

    /* Set data */
    efcm->SetData        ( X    );
    efcm->SetTrueClusters( NULL );
    efcm->SetMaxValue    ( GetMaxValue() );

    /* Run ensemble fuzzy c-means clustering */
    efcm->EFCM           ( logfile );

    /* get result of reclustering */
    int FinalClusters[TotalClusters];
    FC = efcm->GetFinalClusters(FinalClusters);

    /* make cluster prototypes */
    /* Allocate necessary storage to **V (centroids) */
    V = new double * [FC];
    for( i = 0; i < FC; i++ )
        V[i] = new double [D];

    efcm->GetClusters( V );

    /* update cluster label for all data (update **BF) after reclustering */
    for( x = 0, j = 0; j < NrBlocks; j++ )
    {
      for( i = 0; i < FC; i++ )
       for( k = x; k < x + BFC[j]; k++ )
         if( FinalClusters[k] == i )
           for( y = 0; y < BL[j]; y++ )
             if( BF[j][y] == k - x )
                 BF[j][y] = i;
      x += BFC[j];
    }

    /* casting final clusters (**BF) into *F for validation */
    F = new int [N];
    for( k = 0, i = 0; i < NrBlocks; i++ )
      for( j = 0; j < BL[i]; j++ )
           F[k++] = BF[i][j];

    /* casting true clusters (**BTC) into *TC for validation */
    TC = new int [NrBlocks * B];
    for( k = 0, i = 0; i < NrBlocks; i++ )
      for( j = 0; j < BL[i]; j++ )
           TC[k++] = BTC[i][j];

    end = clock();
    printf("\nReclustering of %d blocks with %d samples...",BIdx, TotalClusters);
    printf(" [Done in %d.%d seconds]", (end - start) / CLK_TCK, (end - start) % CLK_TCK);
    printf("\n%d clusters discovered.",FC);

    if( fp != NULL )
    {
        fprintf(fp, " (%d clusters discovered)",FC);
        fprintf(fp, " [Done in %d.%d seconds]", (end - start) / CLK_TCK, (end - start) % CLK_TCK);
    }

    if( Validation )
    {
        double dValidation = ClusterValidation();
        start = clock();
        printf("\nCluster Validity (Jaccard Index): %f ", dValidation);
        end = clock();
        printf(" [Done in %d.%d seconds]", (end - start) / CLK_TCK, (end - start) % CLK_TCK);

        /* log all activities including time consuming */
        if( fp != NULL )
        {
            fprintf(fp, "\nCluster Validity (Jaccard Index): %f ", dValidation);
            fprintf(fp, " [Done in %d.%d seconds]", (end - start) / CLK_TCK, (end - start) % CLK_TCK);
        }
    }    

    fclose(fp);

    delete efcm;

    return 0;
}

/* Cluster Validation =====================================================*/

double StreamEnssembleFuzzyCMeans::ClusterValidation( )
{
    int i, j;

    /* Allocate necessary storage to **Clusters */
    ConfusionMatrix = new long * [K];
    for( i = 0; i < K; i++ ) 
         ConfusionMatrix[i] = new long [FC];

    for( i = 0; i < K; i++ )
      for( j = 0; j < FC; j++ ) 
         ConfusionMatrix[i][j] = 0;

    for( i = 0; i < K; i++ )
      for( j = 0; j < FC; j++ )
      {
           for( long x = 0; x < N; x++ )
             if( TC[x] == i + S )
             {
               for( long y = 0; y < N; y++ )
                 if( F[y] == j && x == y )
                     ConfusionMatrix[i][j]++;
             }
      }

    long double N11, N10, N01, N00;
    /* The category N11 contains the pairs of points that are in the same 
       cluster both in C and in C'. The category N10 contains the pairs of 
       points that are in the same cluster in C but not in C'. The 
       definitions of N01 and N00 are similar. */

    N11 = N10 = N01 = N00 = 0;

    // ----------N11
    for( i = 0; i < K; i++ )
      for( j = 0; j < FC; j++ )
        N11 += ConfusionMatrix[i][j] * ConfusionMatrix[i][j];
    N11 -= N;
    N11 /= 2;

    // ----------N10
    long n = 0;
    for( j = 0; j < FC; j++ )
    {
      n = 0;
      for( i = 0; i < K; i++ )
         n += ConfusionMatrix[i][j];
      N10 += n * n;
    }

    n = 0;
    for( i = 0; i < K; i++ )
      for( j = 0; j < FC; j++ )
         n += ConfusionMatrix[i][j] * ConfusionMatrix[i][j];

    N10 -= n;
    N10 /= 2;

    // ----------N01
    n = 0;
    for( i = 0; i < K; i++ )
    {
      n = 0;
      for( j = 0; j < FC; j++ )
         n += ConfusionMatrix[i][j];
      N01 += n * n;
    }

    n = 0;
    for( i = 0; i < K; i++ )
      for( j = 0; j < FC; j++ )
         n += ConfusionMatrix[i][j] * ConfusionMatrix[i][j];

    N01 -= n;
    N01 /= 2;

    // ----------N00
    long double PP;
    double Jaccard;
    /* PP = N(N - 1)/2 is the total number of point pairs */

    PP = N * ( N - 1 ) / 2;
    N00 = PP - N11 - N10 - N01;

    Jaccard = N11 / ( N11 + N01 + N10 );

    for( i = 0; i < K; i++ )
         delete [] ConfusionMatrix[i];
    delete [] ConfusionMatrix;
    ConfusionMatrix = 0;

    return Jaccard;
}

/* Loading Data ===========================================================*/

void StreamEnssembleFuzzyCMeans::SetData( double **x )
{
        X = x;
}

/*=========================================================================*/

double ** StreamEnssembleFuzzyCMeans::GetData(  )
{
        return X;
}

/*=========================================================================*/

void StreamEnssembleFuzzyCMeans::SetTrueClusters( int * tc )
{
        TC = tc;
}

/*=========================================================================*/

int * StreamEnssembleFuzzyCMeans::GetTrueClusters( )
{
        return TC;
}

/*=========================================================================*/

int * StreamEnssembleFuzzyCMeans::GetTrueClusters( int * tc )
{
        for( long i = 0; i < B; i++ )
            *(tc+i) = *(TC+i);
        return TC;
}

/*=========================================================================*/

void StreamEnssembleFuzzyCMeans::SetMaxValue( double * maxvalue )
{
        MaxValue = maxvalue;
}

/*=========================================================================*/

double * StreamEnssembleFuzzyCMeans::GetMaxValue()
{
        return MaxValue;
}

/* Writing Out Methods ====================================================*/

int StreamEnssembleFuzzyCMeans::WriteMembers( char * FileName )
{
    FILE *fp;
    char  buf[1024];
    long  i;

    clock_t start, end;
    start = clock();
      
    sprintf( buf, "%s.members", FileName );
    printf( "\nLogging %s...", buf );
    if( (fp = fopen( buf, "w" )) == 0 )
    {
        fprintf(stderr,"%s:%d\nError opening %s for mode 'w'... [Fail]\n", __FILE__, __LINE__, buf);
        return 0;
    }

    for( i = 0; i < N; i++ )
         fprintf(fp, "%d\n", F[i]);

    fclose(fp);
    end = clock();
    printf(" [Done in %d.%d seconds]",
    (end - start) / CLK_TCK, (end - start) % CLK_TCK);

    return 0;
}

/*=========================================================================*/

int StreamEnssembleFuzzyCMeans::WriteMRIMembers( char * FileName )
{
    FILE *fp;
    char  buf[1024];
    long  i, k = 0;

    clock_t start, end;
    start = clock();
      
    sprintf( buf, "%s.mri", FileName );
    printf( "\nLogging %s...", buf );
    if( (fp = fopen( buf, "w" )) == 0 )
    {
        fprintf(stderr,"%s:%d\nError opening %s for mode 'w'... [Fail]\n", __FILE__, __LINE__, buf);
        return 0;
    }

    for( i = 0; i < N; i++ )
    {
      fprintf(fp, "%3d", F[i]);
      if( k++ == 255)
      {
          fprintf(fp, "\n");
          k = 0;
      }
    }

    fclose(fp);
    end = clock();
    printf(" [Done in %d.%d seconds]",
    (end - start) / CLK_TCK, (end - start) % CLK_TCK);

    return 0;
}

/*=========================================================================*/

int StreamEnssembleFuzzyCMeans::WriteCentroids( char * FileName )
{
    FILE *fp;
    char  buf[1024];
    int   i, k;
    long  j;

    clock_t start, end;
    start = clock();
      
    sprintf( buf, "%s.centroids", FileName );
    printf( "\nLogging %s...", buf );
    if( (fp = fopen( buf, "w" )) == 0 )
    {
        fprintf(stderr,"%s:%d\nError opening %s for mode 'w'... [Fail]\n", __FILE__, __LINE__, buf);
        return 0;
    }

    for( i = 0; i < FC; i++ )
    {
      for( j = 0; j < D; j++ )
        if( V[i][j] != MAX )
            fprintf(fp, "%f\t", V[i][j]);
      fprintf(fp, "\n");
    }

    fclose(fp);
    end = clock();
    printf(" [Done in %d.%d seconds]",
    (end - start) / CLK_TCK, (end - start) % CLK_TCK);

    return 0;
}

/*=========================================================================*/
