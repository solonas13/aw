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
#include <sys/time.h>
#include <string>
#include "awdefs.h"
#include <stdexcept> // for exceptions
#include <iostream>
#include <cassert>
#include <stack>
#include <utility>
#include <streambuf>
#include <sdsl/suffix_trees.hpp>
#include <sdsl/util.hpp>
#include <sdsl/cst_sct3.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/sdsl_concepts.hpp>
#include <cmath>
#include <string>
using namespace sdsl;
using namespace std;
typedef cst_sct3<> cst_t;
typedef cst_sada<> csts_t;
typedef bp_interval<> node_type;


cst_t cst;
csts_t csts;
INT nnnn;

INT DFSnum2;
INT numk = 0;
INT countleaf=0;
long double numt,numstd = 0;
node_type suffixlinknode;

node_type *DFSunvisited2;
FILE * out_fd;
char * coutfd;

unsigned int compute_aw ( unsigned char * seq, unsigned char * seq_id, struct TSwitch sw )
{
	double start, end;
	INT n = strlen ( ( char * ) seq );
	coutfd=(char*)calloc(n,sizeof(char));
	DFSunvisited2=(node_type*)calloc(n,sizeof(node_type));
	numk = sw.k;
	numt =sw.t;
	nnnn = n;


	if ( ! ( out_fd = fopen ( sw . output_filename, "a") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", sw . output_filename );
		return ( 1 );
	}

    	/* Print the header */
        fprintf ( out_fd, ">%s\n", ( char * ) seq_id );

        fprintf ( out_fd, ".............................................\n");

        fprintf ( out_fd, "k = %ld \n", numk );

        fprintf ( out_fd, "t = %LF \n", numt );

        fprintf ( out_fd, "Avoided Words: \n" );

        /* Compute the Compressed Suffix Tree */
        start = gettime();
        construct_im(cst, (const char *) seq, 1);
        end = gettime();
        fprintf( stderr, " Compressed Suffix Tree construction: %lf secs\n", end - start);

        start = gettime();
        compute_avoidnumk(cst.root(), numk, seq);
        end = gettime();
        fprintf( stderr, " Avoided Words computation: %lf secs\n", end - start);

	if ( fclose ( out_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}
	remove( INPUT_STR );
       	return ( 1 );
       	
       	free(coutfd);
       
       	free(DFSunvisited2);   
}


/* computes the reverse complement of str */
unsigned int RevComStr ( unsigned char * str, unsigned char * str2, INT iLen )
{
   INT i = 0;
   while ( iLen -- )
    {
      switch ( str[iLen] )
       {
         case 'A':
           str2[i++] = 'T';
           break;
         case 'C':
           str2[i++] = 'G';
           break;
         case 'G':
           str2[i++] = 'C';
           break;
         case 'T':
           str2[i++] = 'A';
           break;
         case 'N':
           str2[i++] = 'N';
           break;
         default:
           return ( 0 );
       }
    }
   return ( 1 );
}

double gettime( void )
{
    struct timeval ttime;
    gettimeofday( &ttime , 0 );
    return ttime.tv_sec + ttime.tv_usec * 0.000001;
};
	
	

inline void compute_avoidnumk(const node_type &v, INT nuk, unsigned char * seq){
	
	auto Node = v;
	
	numk = nuk;
	
	DFSnum2 = 0;	
	
	for(auto numdegree = cst.degree(Node); numdegree>0; numdegree--){
		
		   auto NodeChild = cst.select_child(Node, numdegree);
		   
		   DFSnum2++;
		       
		   DFSunvisited2[DFSnum2]=NodeChild;
		   
	}
	
	while(DFSnum2!=0){
		
		auto DFSv2 = DFSunvisited2[DFSnum2];
			
		DFSnum2--;
		
		if(cst.is_leaf(DFSv2)==1){
			
			if(numk<(cst.depth(DFSv2))&&numk>(cst.depth(cst.parent(DFSv2)))){find_avoid(DFSv2, seq);}
			
			}
			
			else{
				
				if(numk<=(cst.depth(DFSv2))&&numk>(cst.depth(cst.parent(DFSv2)))){find_avoid(DFSv2, seq);}
					
					else{
						
						for(auto numdegree = cst.degree(DFSv2); numdegree>0; numdegree--){
		
		           auto NodeChild = cst.select_child(DFSv2, numdegree);
		
		           DFSnum2++;
		       
		           DFSunvisited2[DFSnum2]=NodeChild;}
						
						}
				}	  
}		
}

	
inline void find_avoid(const node_type &v, unsigned char * seq){	
	
	
	long double countnodenow = nodenow(v);
	
	long double countnodesuffix = nodesuffix(v);
	
	long double countnodeprefix = nodeprefix(v);
	
	long double countnodeinfix = nodeinfix(v);
	
	long double max;
	
	long double sign = sqrt((countnodeprefix*countnodesuffix)/countnodeinfix);
	
	numstd = 0.0;
	
	if(sign > 1){max=sign;}else{max=1;}	
		
		long double result = countnodenow-((countnodeprefix*countnodesuffix)/countnodeinfix);
		
		if(result==0){numstd=0;}
			
			else{numstd = result/max;}	
				
		 if(numt >= numstd){	
		 	
		 	INT findposchar=cst.sn(cst.leftmost_leaf(v));
	
	     for(INT f=0; f<numk;f++){coutfd[f]=seq[f+findposchar];}		
	      
	     fprintf ( out_fd, "%s....", (char *) coutfd );
	      
	     fprintf ( out_fd, "std: %LF\n", numstd );	    
	  
	}	
	
	
}

inline INT nodenow(const node_type &v){

  countleaf = 0;
  
  INT herenodenow = cst.size(v);
  
  return herenodenow;
	
}

inline INT nodesuffix(const node_type &v){	
	
	auto Node = cst.sl(v);
	
	countleaf = 0;
	
	INT herenodesuffix;
	
	if(cst.depth(Node)>=(numk-1)&&(numk-1)>(cst.depth(cst.parent(Node)))){
		
		herenodesuffix = cst.size(Node);

		}
		else{ 
			
			do{Node = cst.parent(Node);}while(cst.depth(cst.parent(Node))>=(numk-1));
			herenodesuffix = cst.size(Node);

			}
	
	return herenodesuffix;
		
}

inline INT nodeprefix(const node_type &v){
	
	auto Node = v;
	
	countleaf = 0;
	
	INT herenodeprefix;
	
	if(cst.depth(Node)>=(numk-1)&&(numk-1)>(cst.depth(cst.parent(Node)))){
		
		herenodeprefix = cst.size(Node);
		
		}
		else{ 
						
			herenodeprefix = cst.size(cst.parent(Node));
			}
			
			return herenodeprefix;
	
}

inline INT nodeinfix(const node_type &v){
	
	auto Node = cst.sl(v);
	
	countleaf = 0;
	
	INT herenodeinfix;
	
	
  if(cst.depth(Node)>=(numk-2)&&(numk-2)>(cst.depth(cst.parent(Node)))){
		
		herenodeinfix = cst.size(Node);
		
		}
		else{ 
			
			do{Node = cst.parent(Node);}while(cst.depth(cst.parent(Node))>=(numk-2));
			herenodeinfix = cst.size(Node);
			}
	
	return herenodeinfix;	
}
