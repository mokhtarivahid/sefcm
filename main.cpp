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

#include <iostream>
#include "efcm.h"

#define ON               1
#define OFF              0

/*==========================================================================         */

// Default parameters
static char    FileName[256] = "test.dat";
static int     DataSet       = STATIC;      /* Dataset type(stream/static)           */
static long    BlockSize     = 4*1024;      /* Size of data blocks should be reading */
static double  MinEpsilon    = 0.00001;
static double  MaxEpsilon    = 0.00001;
static double  MinFuzzyExp   = 2.0;
static double  MaxFuzzyExp   = 2.0;
static int     MinC          = 2;           /* Minimum cluster number                */
static int     MaxC          = 2;           /* Maximum cluster number                */
static int     D             = 2;           /* Dimension of data                     */
static long    N             = 2;           /* Data length                           */
static int     C             = 1;           /* Number of clustering ensemble         */
static int     TC            = 2;           /* Number of true clusters in dataset    */
static int     FirstCL       = 0;           /* First cluster label in range dataset  */
static int     ForceC        = 0;           /* Force clustering ensemble to produce 
                                               'ForceC' number clusters              */
static int     Consensus     = COASSOCIATION;/* Consensus strategy                   */
static int     WriteOut      = ON;          /* Write out results                     */
static int     Scan          = ON;          /* Scan data for finding dimension and
                                               length of data                        */
static int     Validation    = OFF;         /* validate cclustering result           */
static char    Method[10]    = "fcm";       /* specify the clustering method that can
                                               be either fcm or efcm                 */
       int     ClusterGeneration = 0;       /* specify how to generate cluster number 
                                               in each run of fcm in efcm, that can be
                                               either 0:random or 1:sequential.      */
       int     OS            = 0;           /* specify the operating system.
                                               0: windows, 1: linux                  */

/*==========================================================================         */

// This function print how usage this program
static void _printUsage(void)
{
       printf("----------------------------------------------------------------------------\n"              );
       printf("-f <filename>   Specifies the input data-file name (default %s).\n"             , FileName   );
       printf("-m <method>     Specifies the clustering method that can be either fcm\n"                    );
       printf("                (fuzzy cmeans) or efcm (ensemble fuzzy cmeans) (default %s).\n" , Method     );
       printf("-s <#>          Specifies the dataset type (default %d).\n"                     , DataSet    );
       printf("                  0: Static Dataset\n"                                                       );
       printf("                  1: Stream Dataset\n"                                                       );
       printf("-b <#>          Specifies the size of data block should be loaded\n"                         );
       printf("                (default %d).\n"                                                , BlockSize  );
       printf("-tc <#>         Specifies the number of true clusters in dataset (default %d).\n", TC        );
       printf("-min_e <#>      Specifies the minimum square error(epsilon) threshold\n"                     );
       printf("                (default %f).\n"                                                , MinEpsilon );
       printf("-max_e <#>      Specifies the maximum square error(epsilon) threshold\n"                     );
       printf("                (default %f).\n"                                                , MaxEpsilon );
       printf("-min_m <#>      Specifies the minimum fuzzification exponent\n"                              );
       printf("                (default %3.3f).\n"                                             , MinFuzzyExp);
       printf("-max_m <#>      Specifies the maximum fuzzification exponent\n"                              );
       printf("                (default %3.3f).\n"                                             , MaxFuzzyExp);
       printf("-min_c <#>      Specifies the minimum number of clusters (default %d).\n"       , MinC       );
       printf("-max_c <#>      Specifies the maximum number of clusters (default %d).\n"       , MaxC       );
       printf("-d <#>          Specifies the dimension of data (default %d).\n"                , D          );
       printf("-n <#>          Specifies the number(length) of data (default %d).\n"           , N          );
       printf("-c <#>          Specifies the number of clusterings(ensembles) (default %d).\n" , C          );
       printf("-force <#>      Force clustering ensemble to produce specified clusters\n"                   );
       printf("                (default %d).\n"                                                , ForceC     );
       printf("-o <off|on>     Write cluster centers and memberships out (default %d).\n"      , WriteOut   );
       printf("-v <off|on>     Validate clustering result base on Jaccard index (default %d).\n", Validation);
       printf("-i <#>          Specifies which integration strategy to use (default %d)\n"     , Consensus  );
       printf("                  0: Relabeling strategy\n"                                                  );
       printf("                  1: Co-Association strategy\n"                                              );
       printf("-scan <off|on>  If you don't know exactly number\\dimension of input data.\n"                );
       printf("                In this situation program scan data and automaticaly\n"                      );
       printf("                finds length(n) and dimension(d) of data (default %d).\n"       , Scan       );
       printf("                Notic: instances should be in rows and dimensions should\n"                  );
       printf("                be in columns.\n"                                                            );
       printf("-os <#>         Specifies the operating system (default %d).\n"                 , OS         );
       printf("                  0: MSWindows OS\n"                                                         );
       printf("                  1: Linux OS\n"                                                             );
       printf("-fcl <#>        Specifies the first cluster label for Jaccard (default %d).\n"  , FirstCL    );
       printf("-gen <#>        Specify how to generate cluster number in each run of fcm\n"                 );
       printf("                in efcm module, that can be either 0:random or 1:sequential\n"               );
       printf("                (default %d).\n"                                                , ClusterGeneration );
       printf("-help|-h        Print how usage this program.\n"                                             );
       printf("----------------------------------------------------------------------------\n"              );
}

/*=========================================================================*/

// This function write arguments into param.ini
static void _writeArgs()
{
       FILE *fp;
       if( (fp = fopen( "sefcm.ini", "w" )) != 0 )
       {
           fprintf(fp, "[FileName]\n%s\n"        , FileName   );
           fprintf(fp, "[ClusteringMethod]\n%s\n", Method     );
           fprintf(fp, "[DataSetType]\n%d\n"     , DataSet    );
           fprintf(fp, "[BlockSize]\n%d\n"       , BlockSize  );
           fprintf(fp, "[TrueClusters]\n%d\n"    , TC         );
           fprintf(fp, "[MinimumEpsilon]\n%f\n"  , MinEpsilon );
           fprintf(fp, "[MaximumEpsilon]\n%f\n"  , MaxEpsilon );
           fprintf(fp, "[MinimumFuzzyExp]\n%f\n" , MinFuzzyExp);
           fprintf(fp, "[MaximumFuzzyExp]\n%f\n" , MaxFuzzyExp);
           fprintf(fp, "[MinimumClusters]\n%d\n" , MinC       );
           fprintf(fp, "[MaximumClusters]\n%d\n" , MaxC       );
           fprintf(fp, "[DataLength]\n%d\n"      , N          );
           fprintf(fp, "[Dimensions]\n%d\n"      , D          );
           fprintf(fp, "[Clusterings]\n%d\n"     , C          );
           fprintf(fp, "[ForceClusters]\n%d\n"   , ForceC     );
           fprintf(fp, "[ConsensusMethod]\n%d\n" , Consensus  );
           fprintf(fp, "[WriteOut]\n%d\n"        , WriteOut   );
           fprintf(fp, "[Validation]\n%d\n"      , Validation );
           fprintf(fp, "[ScanData]\n%d\n"        , Scan       );
       }
       fclose(fp);
}

/*=========================================================================*/

// This function print arguments
static void _printArgs()
{
       printf("\n/*==========================Parameters===========================*/\n");
       printf("[FileName]:\t\t%s\n"       , FileName   );
       printf("[ClusteringMethod]:\t%s\n" , Method     );
       printf("[DataSetType]:\t\t%d\n"    , DataSet    );
       printf("[BlockSize]:\t\t%d\n"      , BlockSize  );
       printf("[TrueClusters]:\t\t%d\n"   , TC         );
       printf("[MinimumEpsilon]:\t%f\n"   , MinEpsilon );
       printf("[MaximumEpsilon]:\t%f\n"   , MaxEpsilon );
       printf("[MinimumFuzzyExp]:\t%f\n"  , MinFuzzyExp);
       printf("[MaximumFuzzyExp]:\t%f\n"  , MaxFuzzyExp);
       printf("[MinimumClusters]:\t%d\n"  , MinC       );
       printf("[MaximumClusters]:\t%d\n"  , MaxC       );
       printf("[DataLength]:\t\t%d\n"     , N          );
       printf("[Dimensions]:\t\t%d\n"     , D          );
       printf("[Clusterings]:\t\t%d\n"    , C          );
       printf("[ForceClusters]:\t%d\n"    , ForceC     );
       printf("[ConsensusMethod]:\t%d\n"  , Consensus  );
       printf("[WriteOut]:\t\t%d\n"       , WriteOut   );
       printf("[Validation]:\t\t%d\n"     , Validation );
       printf("[ScanData]:\t\t%d\n"       , Scan       );
       printf("/*===============================================================*/\n");
}

/*=========================================================================*/

// This function process arguments
static void _processArgs(int argc, char *argv[])
{
   /* HERE on the ones that use the next arg make sure it is there */
   for(int i = 1 ; i < argc ; i++) {
      if(!strcmp(argv[i], "-f")) {
         sscanf(argv[i+1], "%s", FileName);
         /* ignore the next argument */
         i++;
      }else if(!strcmp(argv[i], "-m")) {
         sscanf(argv[i+1], "%s", Method);
         /* ignore the next argument */
         i++;
      } else if(!strcmp(argv[i], "-s")) {
         sscanf(argv[i+1], "%d", &DataSet);
         /* ignore the next argument */
         i++;
      } else if(!strcmp(argv[i], "-b")) {
         sscanf(argv[i+1], "%d", &BlockSize);
         /* ignore the next argument */
         i++;
      } else if(!strcmp(argv[i], "-tc")) {
         sscanf(argv[i+1], "%d", &TC);
         /* ignore the next argument */
         i++;
      } else if(!strcmp(argv[i], "-fcl")) {
         sscanf(argv[i+1], "%d", &FirstCL);
         /* ignore the next argument */
         i++;
      } else if(!strcmp(argv[i], "-gen")) {
         sscanf(argv[i+1], "%d", &ClusterGeneration);
         /* ignore the next argument */
         i++;
      } else if(!strcmp(argv[i], "-min_e")) {
         MinEpsilon = atof( argv[i+1] );
         /* ignore the next argument */
         i++;
      } else if(!strcmp(argv[i], "-max_e")) {
         MaxEpsilon = atof( argv[i+1] );
         /* ignore the next argument */
         i++;
      } else if(!strcmp(argv[i], "-min_m")) {
         MinFuzzyExp = atof( argv[i+1] );
         /* ignore the next argument */
         i++;
      } else if(!strcmp(argv[i], "-max_m")) {
         MaxFuzzyExp = atof( argv[i+1] );
         /* ignore the next argument */
         i++;
      } else if(!strcmp(argv[i], "-min_c")) {
         sscanf(argv[i+1], "%d", &MinC);
         /* ignore the next argument */
         i++;
      } else if(!strcmp(argv[i], "-max_c")) {
         sscanf(argv[i+1], "%d", &MaxC);
         /* ignore the next argument */
         i++;
      } else if(!strcmp(argv[i], "-d")) {
         sscanf(argv[i+1], "%d", &D);
         /* ignore the next argument */
         i++;
      } else if(!strcmp(argv[i], "-n")) {
         sscanf(argv[i+1], "%d", &N);
         /* ignore the next argument */
         i++;
      } else if(!strcmp(argv[i], "-c")) {
         sscanf(argv[i+1], "%d", &C);
         /* ignore the next argument */
         i++;
      } else if(!strcmp(argv[i], "-force")) {
         sscanf(argv[i+1], "%d", &ForceC);
         /* ignore the next argument */
         i++;
      } else if(!strcmp(argv[i], "-o")) {
         if(!strcmp(argv[i+1], "off") || !strcmp(argv[i+1], "0"))
         {
           WriteOut = OFF;
           /* ignore the next argument */
           i++;
         }
         else if(!strcmp(argv[i+1], "on") || !strcmp(argv[i+1], "1"))
         {
           WriteOut = ON;
           /* ignore the next argument */
           i++;
         }
      } else if(!strcmp(argv[i], "-v")) {
         if(!strcmp(argv[i+1], "on") || !strcmp(argv[i+1], "1"))
         {
           Validation = ON;
           /* ignore the next argument */
           i++;
         }
         else if(!strcmp(argv[i+1], "off") || !strcmp(argv[i+1], "0"))
         {
           Validation = OFF;
           /* ignore the next argument */
           i++;
         }
      } else if(!strcmp(argv[i], "-i")) {
         sscanf(argv[i+1], "%d", &Consensus);
         /* ignore the next argument */
         i++;
      } else if(!strcmp(argv[i], "-scan")) {
         if(!strcmp(argv[i+1], "off") || !strcmp(argv[i+1], "0"))
         {
           Scan = OFF;
           /* ignore the next argument */
           i++;
         }
         else if(!strcmp(argv[i+1], "on") || !strcmp(argv[i+1], "1"))
         {
           Scan = ON;
           /* ignore the next argument */
           i++;
         }
      } else if(!strcmp(argv[i], "-os")) {
         sscanf(argv[i+1], "%d", &OS);
         /* ignore the next argument */
         i++;
      } else if(!strcmp(argv[i], "-help")) {
         _printUsage();
         _writeArgs();
         exit(0);
      } else if(!strcmp(argv[i], "-h")) {
         _printUsage();
         _writeArgs();
         exit(0);
      } else {
         printf("Unknown argument: %s.  use -h|-help for help\n", argv[i]);
         exit(0);
      }
   }
}

/*=========================================================================*/

// This function read arguments from param.ini
static void _readArgs()
{
  FILE *fp;
  if( (fp = fopen( "sefcm.ini", "r" )) == 0 )
  {
     _writeArgs();
     return;
  }
  char *str = new char [256];
  do{
      fscanf(fp, "%s" , str );
      
      if(!strcmp(str, "[FileName]"))
         fscanf(fp, "%s" , FileName);
      else if(!strcmp(str, "[ClusteringMethod]"))
         fscanf(fp, "%s" , Method);
      else if(!strcmp(str, "[DataSetType]"))
         fscanf(fp, "%d" , &DataSet);
      else if(!strcmp(str, "[BlockSize]"))
         fscanf(fp, "%d" , &BlockSize);
      else if(!strcmp(str, "[TrueClusters]"))
         fscanf(fp, "%d" , &TC);
      else if(!strcmp(str, "[MinimumEpsilon]" ))
      {
         fscanf(fp, "%s" , str);
         MinEpsilon = atof(str);
      }
      else if(!strcmp(str, "[MaximumEpsilon]" ))
      {
         fscanf(fp, "%s" , str);
         MaxEpsilon = atof(str);
      }
      else if(!strcmp(str, "[MinimumFuzzyExp]"))
      {
         fscanf(fp, "%s" , str);
         MinFuzzyExp = atof(str);
      }
      else if(!strcmp(str, "[MaximumFuzzyExp]"))
      {
         fscanf(fp, "%s" , str);
         MaxFuzzyExp = atof(str);
      }
      else if(!strcmp(str, "[MinimumClusters]"))
         fscanf(fp, "%d" , &MinC);
      else if(!strcmp(str, "[MaximumClusters]"))
         fscanf(fp, "%d" , &MaxC);
      else if(!strcmp(str, "[DataLength]"))
         fscanf(fp, "%d" , &N);
      else if(!strcmp(str, "[Dimensions]"))
         fscanf(fp, "%d" , &D);
      else if(!strcmp(str, "[Clusterings]"))
         fscanf(fp, "%d" , &C);
      else if(!strcmp(str, "[ForceClusters]"))
         fscanf(fp, "%d" , &ForceC);
      else if(!strcmp(str, "[ConsensusMethod]"))
         fscanf(fp, "%d" , &Consensus);
      else if(!strcmp(str, "[WriteOut]"))
         fscanf(fp, "%d" , &WriteOut);
      else if(!strcmp(str, "[Validation]"))
         fscanf(fp, "%d" , &Validation);
      else if(!strcmp(str, "[ScanData]"))
         fscanf(fp, "%d" , &Scan);

  }while( !feof(fp) );
  delete [] str;
  str = 0;
  fclose(fp);
}

/*=========================================================================*/

int main(int argc, char *argv[])
{
  /* Read default parameters from param.ini */
  _readArgs();
  _processArgs(argc, argv);

  if( OS == 0 ) // windows
      system("cls");
  else          // liunx
      system("clear");

  _printArgs();

  /* Create new data object from Load class to loading data */
  Load data;

  /* Scan input file to extract data length(N) and dimension(S) */
  if( Scan == ON )
  {
      if( OS == 0 ) // windows
          system("cls");
      else          // liunx
          system("clear");

      data.ScanData(FileName,N,D);
      printf("\nThe length of data(n)     = %d", N);
      printf("\nThe dimensions of data(d) = %d", D);
      printf("\nDo you accept the length(n) and dimensions(d)? [Y/N] ");
      char ch = getchar();
      if( ch == 'N' || ch == 'n' )
      {
          printf("Please enter data length(n): "   );
          scanf ("%d",&N);
          printf("Please enter data dimension(d): ");
          scanf ("%d",&D);
      }

      /* Read default parameters from param.ini */
      _printArgs();
  }

  if( DataSet == STATIC )
  {
      if( !strcmp(Method, "mri-fcm") )
      {
          /* Load MRImage data */
          //N = (256 * 256);
          data.LoadMRI( FileName, N, D );

          /* Print data loaded by data object */
          //data.PrintData();

          /* Create new fcm object by default parameters */
          FuzzyCMeans fcm( MinEpsilon, MinFuzzyExp, MinC, D, N, TC, Validation, FirstCL );
  
          /* Set data */
          fcm.SetData        ( data.GetData()         );
          fcm.SetTrueClusters( data.GetTrueClusters() );
          fcm.SetMaxValue    ( data.GetMaxValue()     );

          /* Print data loaded by fcm object */
          //fcm.PrintData(1);
  
          clock_t start, end;
          start = clock();
          printf("\nRunning Fuzzy C-Means...");

          /* log all activities including time consuming */
          FILE *fp = NULL;
          char  buf[512];
          sprintf( buf, "%s.log", FileName );
          if( (fp = fopen( buf, "a+" )) == 0 )
              fprintf(stderr,"%s:%d\nError opening %s for mode 'w'... [Fail]\n", __FILE__, __LINE__, buf);
          else
          {
              fprintf(fp, "\n>>------------------------------------------------------------------------------");
              fprintf(fp, "\nFuzzy C-Means on '%s' with %d samples, %d dimension and %d clusters",FileName , N, D, MinC);
              fprintf(fp, "\nRunning Fuzzy C-Means...");
          }

          /* Run fuzzy c-means clustering */
          fcm.FCM(1);

          end = clock();
          printf(" [Done in %d.%d seconds]", (end - start) / CLK_TCK, (end - start) % CLK_TCK);

          /* log elapsed time */
          if( fp != NULL )
              fprintf(fp, " [Done in %d.%d seconds]", (end - start) / CLK_TCK, (end - start) % CLK_TCK);

          if( Validation )
          {
              start = clock();
              printf("\nCluster Validity (Jaccard Index): %f ", fcm.ClusterValidation());
              end   = clock();
              printf(" [Done in %d.%d seconds]", (end - start) / CLK_TCK, (end - start) % CLK_TCK);

              /* log cluster validity */
              if( fp != NULL )
              {
                  fprintf(fp, "\nCluster Validity (Jaccard Index): %f ", fcm.ClusterValidation());
                  fprintf(fp, " [Done in %d.%d seconds]", (end - start) / CLK_TCK, (end - start) % CLK_TCK);
              }
          }

          /* Write out fuzzy c-means clustering results */
          if( WriteOut )
          {
              fcm.WriteCentroids( FileName );
              fcm.WriteUmatrix  ( FileName );
              fcm.WriteMembers  ( FileName );
              //fcm.WriteClusters ( FileName );
              fcm.WriteMRIClusters ( 256, FileName );
          }
      }
      else if( !strcmp(Method, "mri-efcm") )
      {
          /* Load MRImage data */
          //N = (256 * 256);
          data.LoadMRI( FileName, N, D );

          /* Print data loaded by data object */
          //data.PrintData();

          /* Create new efcm object by default parameters */
          EnssembleFuzzyCMeans efcm( MinEpsilon, MaxEpsilon, MinFuzzyExp, MaxFuzzyExp,
                                     MinC, MaxC, N, D, TC, FirstCL, C, Validation, ForceC, Consensus );

          //efcm.PrintParameters();
  
          /* Set data */
          efcm.SetData        ( data.GetData()         );
          efcm.SetTrueClusters( data.GetTrueClusters() );
          efcm.SetMaxValue    ( data.GetMaxValue()     );

          /* Run ensemble fuzzy c-means clustering on MRImage */
          efcm.EFCM( FileName );

          /* Write out final matrix */
          if( WriteOut )
          {
              efcm.WriteMembers    ( FileName );
              efcm.WriteCentroids  ( FileName );
              //efcm.WriteClusters  ( FileName );
              efcm.WriteMRIClusters( 256, FileName );
              efcm.WritePartitions ( FileName );

              /* Write out co-association matrix */
              if( Consensus == COASSOCIATION )
              {
                  //efcm.WriteCoAssocMatrix( FileName );
                  efcm.WriteSLTree       ( FileName );
              }
          }

      }
      else if( !strcmp(Method, "fcm") )
      {
          /* Load data */
          data.LoadData( FileName, N, D );

          /* Print data loaded by data object */
          //data.PrintData();

          /* Create new fcm object by default parameters */
          FuzzyCMeans fcm( MinEpsilon, MinFuzzyExp, MinC, D, N, TC, Validation, FirstCL );
  
          /* Set data */
          fcm.SetData        ( data.GetData()         );
          fcm.SetTrueClusters( data.GetTrueClusters() );
          fcm.SetMaxValue    ( data.GetMaxValue()     );

          /* Print data loaded by fcm object */
          //fcm.PrintData(0);
  
          clock_t start, end;
          start = clock();
          printf("\nRunning Fuzzy C-Means...");

          /* log all activities including time consuming */
          FILE *fp = NULL;
          char  buf[512];
          sprintf( buf, "%s.log", FileName );
          if( (fp = fopen( buf, "a+" )) == 0 )
              fprintf(stderr,"%s:%d\nError opening %s for mode 'w'... [Fail]\n", __FILE__, __LINE__, buf);
          else
          {
              fprintf(fp, "\n>>------------------------------------------------------------------------------");
              fprintf(fp, "\nFuzzy C-Means on '%s' with %d samples, %d dimension and %d clusters",FileName , N, D, MinC);
              fprintf(fp, "\nRunning Fuzzy C-Means...");
          }

          /* Run fuzzy c-means clustering */
          fcm.FCM(1);

          end = clock();
          printf(" [Done in %d.%d seconds]", (end - start) / CLK_TCK, (end - start) % CLK_TCK);

          /* log elapsed time */
          if( fp != NULL )
              fprintf(fp, " [Done in %d.%d seconds]", (end - start) / CLK_TCK, (end - start) % CLK_TCK);

          if( Validation )
          {
              start = clock();
              printf("\nCluster Validity (Jaccard Index): %f ", fcm.ClusterValidation());
              end   = clock();
              printf(" [Done in %d.%d seconds]", (end - start) / CLK_TCK, (end - start) % CLK_TCK);

              /* log cluster validity */
              if( fp != NULL )
              {
                  fprintf(fp, "\nCluster Validity (Jaccard Index): %f ", fcm.ClusterValidation());
                  fprintf(fp, " [Done in %d.%d seconds]", (end - start) / CLK_TCK, (end - start) % CLK_TCK);
              }
          }

          /* Write out fuzzy c-means clustering results */
          if( WriteOut )
          {
              fcm.WriteCentroids( FileName );
              fcm.WriteUmatrix  ( FileName );
              fcm.WriteMembers  ( FileName );
              fcm.WriteClusters ( FileName );
          }
      }
      else
      {
          /* Load data */
          data.LoadData( FileName, N, D );

          /* Print data loaded by data object */
          //data.PrintData();

          /* Create new efcm object by default parameters */
          EnssembleFuzzyCMeans efcm( MinEpsilon, MaxEpsilon, MinFuzzyExp, MaxFuzzyExp,
                                     MinC, MaxC, N, D, TC, FirstCL, C, Validation, ForceC, Consensus );

          //efcm.PrintParameters();
  
          /* Set data */
          efcm.SetData        ( data.GetData()         );
          efcm.SetTrueClusters( data.GetTrueClusters() );
          efcm.SetMaxValue    ( data.GetMaxValue()     );

          /* Run ensemble fuzzy c-means clustering */
          efcm.EFCM( FileName );

          /* Write out final matrix */
          if( WriteOut )
          {
              efcm.WriteMembers   ( FileName );
              efcm.WriteCentroids ( FileName );
              efcm.WriteClusters  ( FileName );
              efcm.WritePartitions( FileName );

              /* Write out co-association matrix */
              if( Consensus == COASSOCIATION )
              {
                  //efcm.WriteCoAssocMatrix( FileName );
                  efcm.WriteSLTree       ( FileName );
              }
          }

      }
  }
  //----------------------------------------------------------------------------
  else //if( DataSet == STREAM )
  {
      if( !strcmp(Method, "mri-efcm") )
      {
          long B = BlockSize;

          /* Create new sefcm object by default parameters */
          StreamEnssembleFuzzyCMeans sefcm( MinEpsilon, MaxEpsilon, MinFuzzyExp, MaxFuzzyExp, MinC, MaxC, 
                                            N, D, BlockSize, TC, FirstCL, C, Validation, ForceC, Consensus, "efcm" );

          /* open input stream file and set file pointer */
          data.OpenStreamData( FileName );

          /* initial sefcm object */
          sefcm.Init(FileName);

          while( data.LoadStreamMRI( FileName, B, D ) )
          {
                 /* Set data */
                 sefcm.SetData        ( data.GetData()         );
                 sefcm.SetTrueClusters( data.GetTrueClusters() );
                 sefcm.SetMaxValue    ( data.GetMaxValue()     );

                 /* Run stream ensemble fuzzy c-means clustering */
                 sefcm.SEFCM( B, FileName );

          };

          if( B > 1 )
          {
                /* Set data */
                sefcm.SetData        ( data.GetData()         );
                sefcm.SetTrueClusters( data.GetTrueClusters() );
                sefcm.SetMaxValue    ( data.GetMaxValue()     );

                /* Run stream ensemble fuzzy c-means clustering */
                sefcm.SEFCM( B, FileName );
          }

          /* Run recluster stream ensemble fuzzy c-means clustering */
          sefcm.ReCluster( FileName );

          /* Write out final matrix */
          if( WriteOut )
          {
                sefcm.WriteMRIMembers( FileName );
                sefcm.WriteCentroids ( FileName );
          }
      }
      else if( !strcmp(Method, "mri-fcm") )
      {
          long B = BlockSize;

          /* Create new sefcm object by default parameters */
          StreamEnssembleFuzzyCMeans sefcm( MinEpsilon, MaxEpsilon, MinFuzzyExp, MaxFuzzyExp, MinC, MaxC, 
                                            N, D, BlockSize, TC, FirstCL, C, Validation, ForceC, Consensus, "fcm" );

          /* open input stream file and set file pointer */
          data.OpenStreamData( FileName );

          /* initial sefcm object */
          sefcm.Init(FileName);

          while( data.LoadStreamMRI( FileName, B, D ) )
          {
                 /* Set data */
                 sefcm.SetData        ( data.GetData()         );
                 sefcm.SetTrueClusters( data.GetTrueClusters() );
                 sefcm.SetMaxValue    ( data.GetMaxValue()     );

                 /* Run stream ensemble fuzzy c-means clustering */
                 sefcm.SEFCM( B, FileName );

          };

          if( B > 1 )
          {
                /* Set data */
                sefcm.SetData        ( data.GetData()         );
                sefcm.SetTrueClusters( data.GetTrueClusters() );
                sefcm.SetMaxValue    ( data.GetMaxValue()     );

                /* Run stream ensemble fuzzy c-means clustering */
                sefcm.SEFCM( B, FileName );
          }

          /* Run recluster stream ensemble fuzzy c-means clustering */
          sefcm.ReCluster( FileName );

          /* Write out final matrix */
          if( WriteOut )
          {
                sefcm.WriteMRIMembers( FileName );
                sefcm.WriteCentroids ( FileName );
          }
      }
      else
      {
          long B = BlockSize;

          /* Create new sefcm object by default parameters */
          StreamEnssembleFuzzyCMeans sefcm( MinEpsilon, MaxEpsilon, MinFuzzyExp, MaxFuzzyExp, MinC, MaxC, 
                                            N, D, BlockSize, TC, FirstCL, C, Validation, ForceC, Consensus, Method );

          /* open input stream file and set file pointer */
          data.OpenStreamData( FileName );

          /* initial sefcm object */
          sefcm.Init(FileName);

          while( data.LoadStreamData( FileName, B, D ) )
          {
                 /* Set data */
                 sefcm.SetData        ( data.GetData()         );
                 sefcm.SetTrueClusters( data.GetTrueClusters() );
                 sefcm.SetMaxValue    ( data.GetMaxValue()     );

                 /* Run stream ensemble fuzzy c-means clustering */
                 sefcm.SEFCM( B, FileName );

          };

          if( B > 1 )
          {
                /* Set data */
                sefcm.SetData        ( data.GetData()         );
                sefcm.SetTrueClusters( data.GetTrueClusters() );
                sefcm.SetMaxValue    ( data.GetMaxValue()     );

                /* Run stream ensemble fuzzy c-means clustering */
                sefcm.SEFCM( B, FileName );
          }

          /* Run recluster stream ensemble fuzzy c-means clustering */
          sefcm.ReCluster( FileName );

          /* Write out final matrix */
          if( WriteOut )
          {
                sefcm.WriteMembers   ( FileName );
                sefcm.WriteCentroids ( FileName );
          }
      }
  }

  printf("\n");

  /* Update new parameters into sefcm.ini */
  _writeArgs();

  return 0;
}

/*=========================================================================*/
