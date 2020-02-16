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

#include "load.h"

#define U(i,j) U[j][i]
#define Co_AssocN        (((N)*(N))-N)/2
#define Co_Assoc(i,j)    Co_Assoc[(((i)*(N))+(j)) - (((i)*(i)+((3*(i))+2))/2)]

/* clustering methods */
#define _FCM             "fcm"
#define _EFCM            "efcm"

/* consensus methods */
#define RELABELING       0
#define COASSOCIATION    1

/* dataset types */
#define STATIC           0
#define STREAM           1

/* cluster type */
#define JOINED           0  // integrate all data with a cluster label in a one file
#define DISJOINED        1  // put data with the same cluster label in separated files

extern  int OS;             /* specify the operating system. (defined in main.cpp)
                               0: windows, 1: linux                        */
extern  int ClusterGeneration;
                            /* specify how to generate cluster number in each run 
                               of fcm, that can be either 0:random or 1:sequential. */

/***************************************************************************/
/*                            Fuzzy C-Means Class                          */
/***************************************************************************/

class FuzzyCMeans
{
  public:
        /* Constructors */
        FuzzyCMeans                    (                          );
        FuzzyCMeans                    ( double Epsilon,
                                         double FuzzyExp,
                                         int    NrClusters,
                                         int    NrDimensions,
                                         long   NrData,
                                         int    NrTrueClusters    = 2,
                                         int    ClusterValidation = 0,
                                         int    FirstClusterLabel = 0,
                                         double *MaxValues        = NULL );
      ~ FuzzyCMeans                    (                          );

        /* Parameters */
        void        ResetParameters    ( double Epsilon,
                                         double FuzzyExp,
                                         int    NrClusters,
                                         int    NrDimensions,
                                         long   NrData,
                                         double *MaxValues = NULL );
        void        GetParameters      ( double &Epsilon,
                                         double &FuzzyExp,
                                         int    &NrClusters,
                                         int    &NrDimensions,
                                         long   &NrData           );
        void        GetParameters      ( char   *FileName,
                                         double &Epsilon,
                                         double &FuzzyExp,
                                         int    &NrClusters,
                                         int    &NrDimensions,
                                         long   &NrData           );

        /* Main functions */
        int         FCM                ( int displayinfo = 1      );
        int         Init               (                          );
        int         UpdateCentroids    (                          );
        double      UpdateUmatrix      (                          );
        double      ClusterValidation  (                          );
        int         GetNrFinalCluster  (                          );
        int         GetFinalClusters   ( int    *  f              );
        int         GetClusters        ( double ** v              );

        /* Utilities */
        int         IsExampleCentroid  ( int k                    );
        double      Distance           ( double * v1, double * v2 );

        /* Public Output Utilities */
        double   ** GetCentroids       (                          );
        double   ** GetUmatrix         (                          );
        int      *  GetMembers         ( int    *  m              );
        int         WriteCentroids     ( char   *  filename       );
        int         WriteUmatrix       ( char   *  filename       );
        int         WriteMembers       ( char   *  filename       );
        int         WriteClusters      ( char   *  filename       );
        int         WriteMRIClusters   ( int col = 256, 
                                         char * filename = NULL   );

        /* Referral **X to external **x data */
        void        SetData            ( double ** x              );
        void        SetTrueClusters    ( int    *  tc             );
        void        SetMaxValue        ( double *  maxvalue       );

        /* Print info into console. */
        void        PrintParameters    ( int display = 1          );
        void        PrintData          ( int display = 1          );

  protected:
        double      Epsilon;                // square error threshold
        double      m;                      // fuzzification exponent
        int         C;                      // number of clusters
        int         D;                      // dimension of input data
        long        N;                      // number of input data
        int         NumberOfIterations;
        double   ** V;                      // cluster prototypes (centers)         V[C][D]
        double   ** U;                      // membership matrix                    U[N][C]
        double   ** X;                      // input data                           X[N][D]
        double   *  MaxValue;               // maximum value of attributes (used in first random init) MaxValue[D]

        /* cluster validation attributes */
        int         K;                      // number of true clusters (labels) ,used in validity
        int         S;                      // first cluster label              ,used in validity
        int      *  TC;                     // true clusters (labels)  TC[N]    ,used in validity
        int         Validation;             // validate cclustering result
        long     ** ConfusionMatrix;        /* Virtually all criteria for comparing clusterings can be computed given
                                               the so-called confusion matrix. Assume that we have two clusterings 
                                               C ={C1,C2, . . . ,CK} and C' = {C'1,C'2 , . . . ,C'K'}. The confusion 
                                               matrix M =(mij) is a K x K' matrix whose ijth element is the number 
                                               of points in the intersection of clusters Ci and C'j , 
                                               i.e., mij = |Ci \ C'j |. */
};

/***************************************************************************/
/*                      Ensemble Fuzzy C-Means Class                       */
/***************************************************************************/

class EnssembleFuzzyCMeans
{
  public:
        /* Constructors */
        EnssembleFuzzyCMeans               (                           );
        EnssembleFuzzyCMeans               ( double MinEpsilon        = 0.0001,
                                             double MaxEpsilon        = 0.0001,
                                             double MinFuzzyexp       = 2.0,
                                             double MaxFuzzyexp       = 3.0,
                                             int    MinCluster        = 2,
                                             int    MaxCluster        = 2,
                                             long   NrData            = 2,
                                             int    NrDimension       = 2,
                                             int    NrTrueClusters    = 2,
                                             int    FirstClusterLabel = 0,
                                             int    NrClusterings     = 2,
                                             int    ClusterValidation = 0,
                                             int    NrForceClusters   = 0,
                                             int    ConsensusMethod   = COASSOCIATION );
      ~ EnssembleFuzzyCMeans               (                           );

        /* Parameters */
        void           ResetParameters     ( double MinEpsilon        = 0.0001,
                                             double MaxEpsilon        = 0.0001,
                                             double MinFuzzyexp       = 2.0,
                                             double MaxFuzzyexp       = 3.0,
                                             int    MinCluster        = 2,
                                             int    MaxCluster        = 2,
                                             long   NrData            = 2,
                                             int    NrDimension       = 2,
                                             int    NrTrueClusters    = 2,
                                             int    FirstClusterLabel = 0,
                                             int    NrClusterings     = 2,
                                             int    ClusterValidation = 0,
                                             int    NrForceClusters   = 0,
                                             int    ConsensusMethod   = COASSOCIATION );
        void           GetParameters       ( double &min_epsilon,
                                             double &max_epsilon,
                                             double &min_fuzzyexp,
                                             double &max_fuzzyexp,
                                             int &min_k, int &max_k,
                                             int &d, long &n, int &c   );

        /* Main functions */
        int            EFCM                ( char * logfile = NULL     );
        int            Init                (                           );
        int            MakeClusters        ( char * logfile = NULL     );
        int         *  UpdateCoAssocUmatrix( int c = 0, char * logfile = NULL );
        double         ClusterValidation   (                           );
        int            GetNrFinalCluster   (                           );
        int            GetFinalClusters    ( int * f                   );
        int            GetClusters         ( double ** v               );

        /* Referral **X to external **x data */
        void           SetData             ( double **x                );
        void           SetTrueClusters     ( int    * tc               );
        void           SetMaxValue         ( double * maxvalue         );
        double      *  GetMaxValue         (                           );

        /* Print info into console. */
        void           PrintParameters     ( int display = 1           );
        int            WriteMembers        ( char * filename           );
        int            WriteCentroids      ( char * filename           );
        int            WriteClusters       ( char * filename, 
                                             int integration = JOINED  );
        int            WriteMRIClusters    ( int col = 256, 
                                             char * filename = NULL    );
        int            WritePartitions     ( char * filename           );
        int            WriteCoAssocMatrix  ( char * filename           );
        int            WriteSLTree         ( char * filename           );

  protected:

        /* public attributes */
        double         Min_Epsilon;      // minimum square error threshold
        double         Max_Epsilon;      // maximum square error threshold
        double         Min_m;            // maximum fuzzification exponent
        double         Max_m;            // maximum fuzzification exponent
        int            Min_K;            // maximum number of clusters
        int            Max_K;            // maximum number of clusters
        int            D;                // dimension of input data
        long           N;                // length of input data
        int            C;                // number of clusterings (ensembles)
        int            FC;               // number of final clusters
        double      ** X;                // input data                    X[N][D]
        int         *  F;                // final partition (cluster name of each element) F[N]
        double      ** V;                // cluster prototypes (centers)  V[FC][D]
        int            Consensus;        // integration strategy (consensus)
        double      *  MaxValue;         // maximum value of attributes   MaxValue[D]
        FuzzyCMeans *  fcm;

        /* relabeling & voting attributes */
        int         ** P;                // partitions                    P[C][N] (used in relabeling & voting matrix)

        /* evidence accumulation attributes */
        int         ** FinalClusters;    // FinalClusters[C][N]
        int            ForceC;           // number of final clusters which efcm forced to produce them
        int         *  Co_Assoc;         // co-association matrix         C[N][N] in fact: (C[(N*N-N)/2])
        long        ** SL[2];            /* single-linkage tree           SL[2][N][C]
                                            first dimension is either '0' (left) or '1' (right) which means '0'
                                            contains left neighbor and '1' contains right neighbor in its cluster.
                                            second dimension is the index of data in **X for left and right neighbor.
                                            third dimension is the level of single-linkage tree.
                                            every SL[][i][j] contains the left and right co-clusters for i'th data
                                            in j'th level of single-linkage tree. */
        long        ** Clusters;         /* store number of clusters and level of each partition. 
                                            first dimension contains two values involve 
                                            cluster-number ([0]) and level of partition ([1]).
                                            second dimension is the maximum partition number ([C+1]),
                                            actually is the number of ensembles */
        long           Clustersidx;      /* index of Clusters[2][C+1] (used in second dimension) */

        /* cluster validation attributes */
        int            K;                // number of true clusters (labels) ,used in validity
        int            S;                // first cluster label              ,used in validity
        int         *  TC;               // true clusters (labels)  TC[N]    ,used in validity
        int            Validation;       // validate cclustering result
        long        ** ConfusionMatrix;  /* Virtually all criteria for comparing clusterings can be computed given
                                            the so-called confusion matrix. Assume that we have two clusterings 
                                            C ={C1,C2, . . . ,CK} and C' = {C'1,C'2 , . . . ,C'K'}. The confusion 
                                            matrix M =(mij) is a K x K' matrix whose ijth element is the number 
                                            of points in the intersection of clusters Ci and C'j , 
                                            i.e., mij = |Ci \ C'j |. */
};

/***************************************************************************/
/*                  Stream Ensemble Fuzzy C-Means Class                    */
/***************************************************************************/

class StreamEnssembleFuzzyCMeans
{
  public:
        /* Constructors */
        StreamEnssembleFuzzyCMeans         (                           );
        StreamEnssembleFuzzyCMeans         ( double MinEpsilon        = 0.0001,
                                             double MaxEpsilon        = 0.0001,
                                             double MinFuzzyexp       = 2.0,
                                             double MaxFuzzyexp       = 3.0,
                                             int    MinCluster        = 2,
                                             int    MaxCluster        = 2,
                                             long   NrData            = 2,
                                             int    NrDimension       = 2,
                                             long   BlockSize         = 2,
                                             int    NrTrueClusters    = 2,
                                             int    FirstClusterLabel = 0,
                                             int    NrClusterings     = 2,
                                             int    ClusterValidation = 0,
                                             int    NrForceClusters   = 0,
                                             int    ConsensusMethod   = COASSOCIATION,
                                             char * BaseClusterMethod = _EFCM );
      ~ StreamEnssembleFuzzyCMeans         (                           );

        /* Parameters */
        void           ResetParameters     ( double MinEpsilon        = 0.0001,
                                             double MaxEpsilon        = 0.0001,
                                             double MinFuzzyexp       = 2.0,
                                             double MaxFuzzyexp       = 3.0,
                                             int    MinCluster        = 2,
                                             int    MaxCluster        = 2,
                                             long   NrData            = 2,
                                             int    NrDimension       = 2,
                                             long   BlockSize         = 1000,
                                             int    NrTrueClusters    = 2,
                                             int    FirstClusterLabel = 0,
                                             int    NrClusterings     = 2,
                                             int    ClusterValidation = 0,
                                             int    NrForceClusters   = 0,
                                             int    ConsensusMethod   = COASSOCIATION,
                                             char * BaseClusterMethod = _EFCM );

        /* Main functions */
        int            Init                ( char * logfile = NULL     );
        int            SEFCM               ( long block = 1000, 
                                             char * logfile = NULL     );
        int            ReCluster           ( char * logfile = NULL     );
        double         ClusterValidation   (                           );

        /* Referral **X to external **x data */
        void           SetData             ( double **x                );
        double    **   GetData             (                           );
        void           SetTrueClusters     ( int    * tc               );
        int       *    GetTrueClusters     (                           );
        int       *    GetTrueClusters     ( int    * tc               );
        void           SetMaxValue         ( double * maxvalue         );
        double    *    GetMaxValue         (                           );

        /* Print info into console. */
        void           PrintParameters     ( int  display = 1          );
        int            WriteMembers        ( char * filename           );
        int            WriteMRIMembers     ( char * filename           );
        int            WriteCentroids      ( char * filename           );
        int            WriteClusters       ( char * filename           );

  protected:

        /* public attributes */
        double         Min_Epsilon;      // minimum square error threshold
        double         Max_Epsilon;      // maximum square error threshold
        double         Min_m;            // maximum fuzzification exponent
        double         Max_m;            // maximum fuzzification exponent
        int            Min_K;            // maximum number of clusters
        int            Max_K;            // maximum number of clusters
        int            D;                // dimension of input data
        long           N;                // length of input data
        long           B;                // size of blocks should be loaded from input data iteratively
        int            C;                // number of clusterings (ensembles)
        int            FC;               // number of final clusters
        double    **   X;                // input data                    X[N][D]
        int       *    F;                // final partition (cluster name of each element) F[N]
        double    **   V;                // cluster prototypes (centers)  V[FC][D]
        int            Consensus;        // integration strategy (consensus)
        double    *    MaxValue;         // maximum value of attributes   MaxValue[D]
        char      *    BaseMethod;       // base clustering method which can be either "efcm" or "fcm"

        /* evidence accumulation attributes */
        int            ForceC;           // number of final clusters which efcm forced to produce them

        /* stream attributes at every block */
        int            BIdx;             // blocks indexes that proceed to 'Blocks'
        int            NrBlocks;         // block numbers that equales to (N / B)
        int       *    BFC;              // number of final clusters at every block                       FC[NrBlocks]
        int       **   BF;               // final partition at every block (cluster name of each element) F[NrBlocks][B]
        double    ***  BV;               // cluster prototypes (centers) at every block                   V[NrBlocks][BFC][D]
        long      *    BL;               // data length for each efcm clustering at every block           L[NrBlocks]
        int       **   BTC;              // true clusters (labels) at every block                         TC[NrBlocks][B]

        /* cluster validation attributes */
        int            K;                // number of true clusters (labels) ,used in validity
        int            S;                // first cluster label              ,used in validity
        int         *  TC;               // true clusters (labels)  TC[N]    ,used in validity
        int            Validation;       // validate cclustering result
        long        ** ConfusionMatrix;  /* Virtually all criteria for comparing clusterings can be computed given
                                            the so-called confusion matrix. Assume that we have two clusterings 
                                            C ={C1,C2, . . . ,CK} and C' = {C'1,C'2 , . . . ,C'K'}. The confusion 
                                            matrix M =(mij) is a K x K' matrix whose ijth element is the number 
                                            of points in the intersection of clusters Ci and C'j , 
                                            i.e., mij = |Ci ^ C'j|. */
};
