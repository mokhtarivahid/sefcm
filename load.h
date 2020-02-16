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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#ifndef CLK_TCK
#define CLK_TCK           1000000
#endif
#define MIN              -2147483648.0
#define MAX               2147483647.0
#define MININT           -32768
#define MAXINT            32768
#define UNKNOWNINT       -32768
#define NEWLINE           '\n'
#define SPACE             ' '
#define TAB               '\t'
#define COMMA             ','
#define STREND            '\0'
#define IsDelimiter(X)    ( (X) == TAB || (X) == SPACE || (X) == COMMA  )
#define IsEndOfString(X)  ( (X) == NEWLINE || (X) == STREND )
#define IsEndOfWord(X)    ( (X) == TAB || (X) == SPACE   || (X) == COMMA || (X) == NEWLINE || (X) == STREND )
#define IsComment(X)      ( (X) == '%' || (X) == '/' || (X) == '{' || (X) == '}' || (X) == '[' || (X) == ']' )


/***************************************************************************/
/*                            Load Class                                   */
/***************************************************************************/

class Load
{
  public:
        /* Constructors */
        Load                      (                                 );
        Load                      ( char * filename                 );
        Load                      ( char * filename, int d          );
        Load                      ( char * filename, long n, int d  );
      ~ Load                      (                                 );

        /* Main methods */
        double   ** LoadData      (                                 );
        double   ** LoadData      ( char * filename, long n, int d  );
        double   ** LoadMRI       ( char * filename, long &n,int d  );
        double   ** LoadMRI2      ( char * filename, long &n,int d  );

        /* just open the stream file */
        int         OpenStreamData( char * filename                 );

        /* load data from the stream file */
        int         LoadStreamData( char * filename, long &n, int d );
        int         LoadStreamMRI ( char * filename, long &n, int d );
        int         ScanData      ( char * filename, long &n, int &d);
        double   ** GetData       (                                 );
        double   ** GetData       ( long &n, int &d                 );
        double   *  GetMaxValue   (                                 );
        int      *  GetTrueClusters(                                );
        void        PrintData     ( int display = 1                 );

  protected:
        int         D;         /* dimension of input data                    */
        long        N;         /* length of input data                       */
        double   ** X;         /* input data: X[N][D]                        */
        int      *  TC;        /* true cluster label for each entry: TC[N]   */
        char     *  FileName;  /* input file name                            */
        double   *  MaxValue;  /* maximum value of attributes (used in first 
                                  random init in fcm class) MaxValue[D]      */
        FILE     *  sfp;       /* stream file pointer                        */
        int         CurBlock;  /* number of current block in the stream file */
};
