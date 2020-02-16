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

/***************************************************************************/
/*                            Load Class                                   */
/***************************************************************************/

/* Constructors ===========================================================*/

Load::Load()
{
        D                       = 2;
        N                       = 2;
        X                       = NULL;
        TC                      = NULL;
        FileName                = new char [512];
        MaxValue                = NULL;
        CurBlock                = 0;
}

Load::Load( char * filename )
{
        FileName = new char [512];
        strcpy( FileName, filename );

        CurBlock = 0;

        /* scan data for reading length(n) and dimension(d) of data */ 
        ScanData( FileName, N, D );

        /* Allocate necessary storage */
        X = new double * [N];
        for( long i = 0; i < N; i++ ) 
             X[i] = new double [D];

        TC = new int [N];

        MaxValue = new double [D];
        for( int i = 0; i < D; i++ )
             MaxValue[i] = MIN;

        /* load data into X[N][D] */
        LoadData();
}

Load::Load( char * filename, int d )
{
        FileName = new char [512];
        strcpy( FileName, filename );

        CurBlock = 0;

        /* scan data for reading length(n) and dimension(d) of data */ 
        ScanData( FileName, N, D );

        D        = d;

        /* Allocate necessary storage */
        X = new double * [N];
        for( long i = 0; i < N; i++ ) 
             X[i] = new double [D];

        TC = new int [N];

        MaxValue = new double [D];
        for( int i = 0; i < D; i++ )
             MaxValue[i] = MIN;

        /* load data into X[N][D] */
        LoadData();
}

Load::Load( char * filename, long n, int d )
{
        FileName = new char [512];
        strcpy( FileName, filename );

        D        = d;
        N        = n;
        CurBlock = 0;

        /* Allocate necessary storage */
        X = new double * [N];
        for( long i = 0; i < N; i++ ) 
             X[i] = new double [D];

        TC = new int [N];

        MaxValue = new double [D];
        for( int i = 0; i < D; i++ )
             MaxValue[i] = MIN;

        /* load data into X[N][D] */
        LoadData();
}

Load::~ Load()
{
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
                TC = NULL;
        }
        if( MaxValue != NULL )
        {
            delete [] MaxValue;
            MaxValue = NULL;
        }
        if( FileName != NULL )
        {
            delete [] FileName;
            FileName = NULL;
        }
}

/* Loading Data ===========================================================*/

double ** Load::LoadData()
{
        clock_t start, end;
        start = clock();

        if( X == NULL )
        {
            fprintf(stderr,"%s:%d\nNo storage allocated to **X !\n", __FILE__, __LINE__);
            exit(1);
        }
        
        FILE *fp;
        fprintf(stderr, "\nLoading data from %s...", FileName);
        if( (fp = fopen( FileName, "r" )) == 0 )
        {
            fprintf(stderr,"%s:%d\nError opening %s for mode 'r'... [Fail]\n", __FILE__, __LINE__, FileName);
            exit(1);
        }

        char *str = new char [1024];
        char *num = new char [50];
        long  row = -1;
        do{
           if( fgets( str, 1024, fp ) != NULL )
           {
               if( row < N - 1 ) row++;
               int i = 0, col = 0;
               while( ! IsEndOfString(str[i]) )
               {
                      while( IsDelimiter(str[i]) ) i++;
                      if( IsComment(str[i]) ) { row--; break; }

                      int j = 0;
                      while( ! IsEndOfWord(str[i]) )
                             num[j++] = str[i++];
                      num[j] = STREND;
                      double n = atof( num );
                      TC[row] = (int)n;
                      if( col < D && MaxValue[col] < n ) MaxValue[col] = n;
                      if( col < D ) X[row][col++] = n;
               }
           }
        }while( !feof(fp) );

        delete [] str;
        delete [] num;
        fclose( fp );

        end = clock();
        printf(" [Done in %d.%d seconds]",
        (end - start) / CLK_TCK, (end - start) % CLK_TCK);

        return X;
}

/*=========================================================================*/

double ** Load::LoadData( char * filename, long n, int d )
{
        strcpy( FileName, filename );

        D        = d;
        N        = n;

        /* Allocate necessary storage */
        if( X != NULL )
        {
                for( int i = 0; i < D; i++ )
                        delete [] X[i];
                delete [] X;
                X = NULL;
        }

        X = new double * [N];
        for( long i = 0; i < N; i++ ) 
             X[i] = new double [D];

        if( TC != NULL )
        {
                delete [] TC;
                TC = NULL;
        }
        TC = new int [N];

        MaxValue = new double [D];
        for( int i = 0; i < D; i++ )
             MaxValue[i] = MIN;

        /* load data into X[N][D] */
        return LoadData();
}

/* MRI data ===============================================================*/

double ** Load::LoadMRI( char * filename, long &n, int d )
{
        clock_t start, end;
        start = clock();

        strcpy( FileName, filename );

        D        = d;
        N        = n;

        /* Allocate necessary storage */
        if( X != NULL )
        {
                for( int i = 0; i < D; i++ )
                        delete [] X[i];
                delete [] X;
                X = NULL;
        }

        X = new double * [N];
        for( long i = 0; i < N; i++ ) 
             X[i] = new double [D];

        MaxValue = new double [D];
        for( int i = 0; i < D; i++ )
             MaxValue[i] = 256;

        FILE *fp;
        fprintf(stderr, "\nLoading MRI image from %s...", FileName);
        if( (fp = fopen( FileName, "r" )) == 0 )
        {
            fprintf(stderr,"%s:%d\nError opening %s for mode 'r'... [Fail]\n", __FILE__, __LINE__, FileName);
            exit(1);
        }

        char *str = new char [1024];
        char *num = new char [50];
        long  row = 0;
        do{
           if( fgets( str, 1024, fp ) != NULL )
           {
               int i = 0;
               while( ! IsEndOfString(str[i]) )
               {
                      while( IsDelimiter(str[i]) ) i++;
                      if( IsComment(str[i]) ) { row--; break; }

                      int j = 0;
                      while( ! IsEndOfWord(str[i]) )
                             num[j++] = str[i++];
                      num[j] = STREND;
                      if( row < N ) X[row++][0] = atoi( num );
               }
           }
        }while( !feof(fp) );

        delete [] str;
        delete [] num;
        fclose( fp );

        end = clock();
        printf(" [Done in %d.%d seconds]",
        (end - start) / CLK_TCK, (end - start) % CLK_TCK);
        
        return X;
}

/* ========================================================================*/

double ** Load::LoadMRI2( char * filename, long &n, int d )
{
        clock_t start, end;
        start = clock();

        strcpy( FileName, filename );

        D        = d;
        N        = n;

        FILE *fp;
        fprintf(stderr, "\nLoading MRI image from %s...", FileName);
        if( (fp = fopen( FileName, "r" )) == 0 )
        {
            fprintf(stderr,"%s:%d\nError opening %s for mode 'r'... [Fail]\n", __FILE__, __LINE__, FileName);
            exit(1);
        }

        /* Allocate necessary storage */

        if( X != NULL )
        {
                for( int i = 0; i < D; i++ )
                        delete [] X[i];
                delete [] X;
                X = NULL;
        }

        X = new double * [N];
        for( long i = 0; i < N; i++ ) 
             X[i] = new double [D];

        MaxValue = new double [D];
        for( int i = 0; i < D; i++ )
             MaxValue[i] = 256;

        long i,j;
        unsigned short int *buf;

        buf = new unsigned short int [N];

        for( i = 0; i < D; i++ )
        {
            fread( buf, 2, N, fp );
            for( j = 0; j < N; j++ )
                X[j][i] = buf[j];
        }

        fclose(fp);

        end = clock();
        printf(" [Done in %d.%d seconds]",
        (end - start) / CLK_TCK, (end - start) % CLK_TCK);

        n = N;

        return X;
}

/* Stream data loading ====================================================*/

int Load::OpenStreamData( char * filename )
{
        strcpy( FileName, filename );

        clock_t start, end;
        start = clock();

        fprintf(stderr, "\nLoading stream data from %s...", FileName);
        if( (sfp = fopen( FileName, "r" )) == 0 )
        {
            fprintf(stderr,"%s:%d\nError opening %s for mode 'r'... [Fail]\n", __FILE__, __LINE__, FileName);
            exit(1);
        }

        end = clock();
        printf(" [Done in %d.%d seconds]",
        (end - start) / CLK_TCK, (end - start) % CLK_TCK);
        return 0;
}

/*=========================================================================*/

/* load data from stream data in size of n, whenever file pointer reach to
   end of file this method return 0, otherwise whenever read data in length
   of n return 1.
   this method should be used in a loop, while it return 1, read all stream 
   data sequentially in blocks of length of n.
*/
int Load::LoadStreamData( char * filename, long &n, int d )
{

        D        = d;
        N        = n;

        /* Allocate necessary storage */
        if( X != NULL )
        {
                for( int i = 0; i < D; i++ )
                        delete [] X[i];
                delete [] X;
                X = NULL;
        }
        X = new double * [N];
        for( long i = 0; i < N; i++ ) 
             X[i] = new double [D];

        if( TC != NULL )
        {
                delete [] TC;
                TC = NULL;
        }
        TC = new int [N];

        MaxValue = new double [D];
        for( int i = 0; i < D; i++ )
             MaxValue[i] = MIN;

        /* load data into X[N][D] */
        clock_t start, end;
        start = clock();

        printf("\nLoading block %d from %s...", CurBlock++, FileName);

        if( X == NULL )
        {
            fprintf(stderr,"%s:%d\nNo storage allocated to **X !\n", __FILE__, __LINE__);
            exit(1);
        }
        
        char *str = new char [1024];
        char *num = new char [50];
        long  row = -1;
        do{
           if( row == N - 1 ) 
           {
               end = clock();
               printf(" [Done in %d.%d seconds (%d samples)]",
               (end - start) / CLK_TCK, (end - start) % CLK_TCK, N);
               return 1;
           }
           if( fgets( str, 1024, sfp ) != NULL )
           {
               if( row < N - 1 ) row++;
               int i = 0, col = 0;
               while( ! IsEndOfString(str[i]) )
               {
                      while( IsDelimiter(str[i]) ) i++;
                      if( IsComment(str[i]) ) { row--; break; }

                      int j = 0;
                      while( ! IsEndOfWord(str[i]) )
                             num[j++] = str[i++];
                      num[j] = STREND;
                      double d = atof( num );
                      TC[row] = (int)d;
                      if( col < D && MaxValue[col] < d ) MaxValue[col] = d;
                      if( col < D ) X[row][col++] = d;
               }
           }
        }while( !feof(sfp) );

        delete [] str;
        delete [] num;
        fclose( sfp );

        /* delete further memory */
        for( long i = row + 1; i < N; i++ )
             delete [] X[i];
        n = N = row + 1;

        end = clock();
        printf(" [Done in %d.%d seconds (%d samples)(EOF)]",
        (end - start) / CLK_TCK, (end - start) % CLK_TCK, N);

        return 0;
}

/*=========================================================================*/

/* load data from stream data in size of n, whenever file pointer reach to
   end of file this method return 0, otherwise whenever read data in length
   of n return 1.
   this method should be used in a loop, while it return 1, read all stream 
   data sequentially in blocks of length of n.
*/
int Load::LoadStreamMRI( char * filename, long &n, int d )
{

        D        = d;
        N        = n;

        /* Allocate necessary storage */
        if( X != NULL )
        {
                for( int i = 0; i < D; i++ )
                        delete [] X[i];
                delete [] X;
                X = NULL;
        }
        X = new double * [N];
        for( long i = 0; i < N; i++ ) 
             X[i] = new double [D];

        if( TC != NULL )
        {
            delete [] TC;
            TC = NULL;
        }
        TC = new int [N];

        MaxValue = new double [D];
        for( int i = 0; i < D; i++ )
             MaxValue[i] = 256;

        /* load data into X[N][D] */
        clock_t start, end;
        start = clock();

        printf("\nLoading block %d from %s...", CurBlock++, FileName);

        if( X == NULL )
        {
            fprintf(stderr,"%s:%d\nNo storage allocated to **X !\n", __FILE__, __LINE__);
            exit(1);
        }

        char *str = new char [1024];
        char *num = new char [50];
        long  row = 0;
        do{
           if( fgets( str, 1024, sfp ) != NULL )
           {
               int i = 0;
               while( ! IsEndOfString(str[i]) )
               {
               if( row == N ) 
               {
                   end = clock();
                   printf(" [Done in %d.%d seconds (%d samples)]",
                   (end - start) / CLK_TCK, (end - start) % CLK_TCK, N);
                   return 1;
               }
                      while( IsDelimiter(str[i]) ) i++;
                      if( IsComment(str[i]) ) { row--; break; }

                      int j = 0;
                      while( ! IsEndOfWord(str[i]) )
                             num[j++] = str[i++];
                      num[j] = STREND;
                      int d = atoi( num );
                      TC[row] = d;
                      if( row < N ) X[row++][0] = d;
               }
           }
        }while( !feof(sfp) );

        delete [] str;
        delete [] num;
        fclose( sfp );

        n = N = row;

        end = clock();
        printf(" [Done in %d.%d seconds (%d samples)(EOF)]",
        (end - start) / CLK_TCK, (end - start) % CLK_TCK, N);

        return 0;
}

/*=========================================================================*/

int Load::ScanData( char *filename, long &n, int &d )
{
        clock_t start, end;
        start = clock();

        strcpy( FileName, filename );
        FILE *fp;
        fprintf(stderr, "\nScanning %s for data numbers(n) and dimensions(d)...", FileName);
        if( (fp = fopen( FileName, "r" )) == 0 )
        {
            fprintf(stderr,"%s:%d\nError opening %s for mode 'r'... [Fail]\n", __FILE__, __LINE__, FileName);
            exit(1);
        }

        char *str = new char [1024];
        char *num = new char [50];
        long  row = -1, col;
        do{
           if( fgets( str, 1024, fp ) != NULL )
           {
               row++;
               int i = 0;
               col = 0;
               while( ! IsEndOfString(str[i]) )
               {
                      while( IsDelimiter(str[i]) ) i++;
                      if( IsComment(str[i]) ) { row--; break; }

                      int j = 0;
                      while( ! IsEndOfWord(str[i]) )
                             num[j++] = str[i++];
                      num[j] = STREND;
                      col++;
                }
           }
        }while( !feof(fp) );

        N = row + 1;
        D = col;

        n = N;
        d = D;

        delete [] str;
        delete [] num;
        fclose( fp );

        end = clock();
        printf(" [Done in %d.%d seconds]",
        (end - start) / CLK_TCK, (end - start) % CLK_TCK);

        return 1;
}

/*=========================================================================*/

double ** Load::GetData()
{
        return X;
}

/*=========================================================================*/

double ** Load::GetData( long &n, int &d )
{
        n = N;
        d = D;
        return X;
}

/*=========================================================================*/

double * Load::GetMaxValue()
{
        return MaxValue;
}

/*=========================================================================*/

int * Load::GetTrueClusters()
{
        return TC;
}

/*=========================================================================*/

void Load::PrintData( int display )
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
