/**
    AW: Avoided Words
    Copyright (C) 2016 Jia Gao and Solon P. Pissis. 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include "awdefs.h"

static struct option long_options[] =
 {
   { "alphabet",                required_argument, NULL, 'a' },
   { "input-file",              required_argument, NULL, 'i' },
   { "output-file",             required_argument, NULL, 'o' },
   { "length",                  required_argument, NULL, 'k' },
   { "threshold",               required_argument, NULL, 't' },
   { "absent",                  required_argument, NULL, 'A' },
   { "reverse",                 required_argument, NULL, 'r' },
   { "common",                  required_argument, NULL, 'c' },
   { "help",                    no_argument,       NULL, 'h' },
   { NULL,                      0,                 NULL, 0   }
 };


/* 
Decode the input switches 
*/
int decode_switches ( int argc, char * argv [], struct TSwitch * sw )
 {
   int          oi;
   int          opt;
   double       val;
   char       * ep;
   int          args;

   /* initialisation */
   sw -> alphabet                       = NULL;
   sw -> input_filename                 = NULL;
   sw -> output_filename                = NULL;
   sw -> k                              = 0;
   sw -> K                              = 0;
   sw -> t                              = 0;
   sw -> r                              = 0;
   sw -> A                              = 0;
   sw -> c                              = 0;
   args = 0;

   while ( ( opt = getopt_long ( argc, argv, "a:i:o:k:t:r:A:c:h", long_options, &oi ) ) != - 1 )
    {
      switch ( opt )
       {
         case 'a':
           sw -> alphabet = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> alphabet, optarg );
           args ++;
           break;

         case 'i':
           sw -> input_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> input_filename, optarg );
           args ++;
           break;

         case 'o':
           sw -> output_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> output_filename, optarg );
           args ++;
           break;

         case 'k':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> k = val;
           break;

         case 'A':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> A = val;
           break;

         case 't':
           val = strtod ( optarg, &ep );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> t = val;
           args ++;
           break;

         case 'r':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> r = val;
           break;

         case 'c':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> c = val;
           break;

         case 'h':
           return ( 0 );
       }
    }

   if ( args < 4 )
     {
       usage ();
       exit ( 1 );
     }
   else
     return ( optind );
 }


/* 
Usage of the tool 
*/
void usage ( void )
 {
   fprintf ( stdout, " Usage: aw <options>\n" );
   fprintf ( stdout, " Standard (Mandatory):\n" );
   fprintf ( stdout, "  -a, --alphabet            <str>     `DNA' for nucleotide  sequences or `PROT'\n"
                     "                                      for protein  sequences. \n" );
   fprintf ( stdout, "  -i, --input-file          <str>     (Multi)FASTA input filename.\n" );
   fprintf ( stdout, "  -o, --output-file         <str>     Output filename.\n" );
   fprintf ( stdout, "  -t, --threshold           <dbl>     The threshold for AWs.\n");
   fprintf ( stdout, " Optional:\n" );
   fprintf ( stdout, "  -k, --length              <int>     Fixed length for AWs (default: no fixed).\n");
   fprintf ( stdout, "  -A, --absent              <int>     `1' to check also for absent avoided words\n"
                     "                                      or `0' otherwise (default: 0).\n"
                     "                                      This option cannot be used with `-c 1'.\n" );
   fprintf ( stdout, "  -c, --common              <int>     `1' to check for common words instead of\n"
                     "                                      avoided or `0' otherwise (default: 0).\n"
                     "                                      This option can only be used with `-k <int>'.\n" );
   fprintf ( stdout, "  -r, --reverse             <int>     `1' to check for the reverse complement or\n"
                     "                                      `0' otherwise (default: 0).\n" );
 }

