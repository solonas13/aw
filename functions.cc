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
#include <sdsl/bit_vectors.hpp>					  // include header for bit vectors
#include "stack.h"

using namespace sdsl;
using namespace std;
typedef cst_sct3<> cst_t;
typedef cst_sada<> csts_t;
typedef bp_interval<> node_type;


cst_t cst;
csts_t csts;
INT nnnn;
INT DFSnum1;
INT DFSnum2;
unsigned int numk = 0;
INT countleaf=0;
long double numt,numstd = 0;
long double numt1,numstd1 = 0;
node_type suffixlinknode;

node_type *DFSunvisited1;
node_type *DFSunvisited2;
node_type *DFSunvisitedall;

FILE * out_fd;
char * coutfd;
char * coutfd1;
char * coutfdall;

unsigned int compute_aw ( unsigned char * seq, unsigned char * seq_id, struct TSwitch sw )
{
	double start, end;
	INT n = strlen ( ( char * ) seq );
	numk = sw.k;
	
	numt1 = sw.t;
	numt  = sw.t;
	nnnn = n;

	if ( ! ( out_fd = fopen ( sw . output_filename, "a") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", sw . output_filename );
		return ( 1 );
	}

	/* Compute the Compressed Suffix Tree */
	start = gettime();
	construct_im(cst, (const char *) seq, 1);
	end = gettime();
	fprintf( stderr, " Compressed Suffix Tree construction: %lf secs\n", end - start);
	
	/* Compute all of AWs */
	if ( sw . k == 0 )
	{
		fprintf ( out_fd, "Seq : %s\n", ( char * ) seq_id );
		fprintf ( out_fd, "t = %LF \n", numt );
		fprintf ( out_fd, ".............................................\n");
		fprintf ( out_fd, "Occuring Avoided Words: \n" );
		coutfdall=(char*)calloc(n,sizeof(char));
		DFSunvisitedall=(node_type*)calloc(n,sizeof(node_type));
		compute_all_of_occuring_avoided_words(cst.root(), seq);	
		free(coutfdall);
		free(DFSunvisitedall); 
		end = gettime();
		fprintf( stderr, " Occuring Avoided Words computation: %lf secs\n", end - start);
	}
	else /* Compute words of fixed length k*/
	{
		fprintf ( out_fd, "Seq : %s\n", ( char * ) seq_id );
		fprintf ( out_fd, "k = %d \n", numk );
		fprintf ( out_fd, "t = %LF \n", numt );

		if ( sw . c == 0 )
		{
			fprintf ( out_fd, ".............................................\n");
			fprintf ( out_fd, "Occuring Avoided Words: \n" );
			start = gettime();
			coutfd=(char*)calloc(n,sizeof(char));
			DFSunvisited2=(node_type*)calloc(n,sizeof(node_type));
			compute_avoidnumk(cst.root(), numk, seq);
			free(coutfd);
			free(DFSunvisited2); 
			end = gettime();
			fprintf( stderr, " Occuring Avoided Words computation: %lf secs\n", end - start);
		}
		else
		{
			fprintf ( out_fd, ".............................................\n");
			fprintf ( out_fd, "Common Words: \n" );
			start = gettime();
			coutfd1=(char*)calloc(n,sizeof(char));
			DFSunvisited1=(node_type*)calloc(n,sizeof(node_type));
			compute_frequencynumk(cst.root(), numk, seq);
			free(coutfd1);
			free(DFSunvisited1); 
			end = gettime();
			fprintf( stderr, " Common Words computation: %lf secs\n", end - start);
		}
	}
	
        if( sw . c == 0 ) 
        { 

		if ( sw . A )
		{
			start = gettime();
			TMaw * Occ = NULL;
			INT NOcc = 0;
			std::reverse(seq, seq + n);
			if ( sw . k == 0 )
			{
				sw . k = 2;
				sw . K = n;
			}
			else
			{
				sw . K = sw . k;
			}
			compute_maw ( seq, seq_id, sw, &Occ, &NOcc );
			std::reverse(seq, seq + n);
			fprintf ( out_fd, ".............................................\n");
			fprintf ( out_fd, "Absent Avoided Words: \n" );

			for ( INT j = 0; j < NOcc; j++ )
			{
			       char * maw = ( char * ) calloc( ( Occ[j] . size + 2 ) , sizeof( char ) );

			       Occ[j] . pos = n - 1 - ( Occ[j] . pos + Occ[j] . size - 1 );
			      
			       auto candidatenode = cst.root();
			       
			       auto candidateinfix = cst.root();
			       
			       auto candidatesuffix = cst.root();
			       
			       long double leavesprefix = 0;
			       
			       long double leavesinfix = 0;
			       
			       long double leavessuffix = 0;
			       
			       /*find prefix node*/
			       
			       if(cst.is_leaf(cst.child(cst.root(),seq[Occ[j]. pos]))==1){
				
				candidatenode = cst.child(cst.root(),seq[Occ[j].pos]);
			
			       }
				
			       else{
			       
			       if(cst.depth(cst.child(cst.root(),seq[Occ[j].pos]))<Occ[j].size){	
			       
			       auto newnode = cst.root();
				
			       auto new2node = cst.root();
			       
			       do{
								
					new2node = cst.child(newnode,(seq[Occ[j].pos+cst.depth(newnode)]));
					
					newnode = new2node;
					
					}
			       while(cst.depth(newnode)<Occ[j].size);
			       
			       if(cst.is_leaf(newnode)==1){
				
				if(cst.depth(newnode)>(Occ[j].size)){
					
					candidatenode = newnode;
					
					leavesprefix = cst.size(candidatenode);

					}
				
				}
				
				else{
					
					candidatenode = newnode;
					
					leavesprefix = cst.size(candidatenode);
					
					}
			       
			       }
			       else{
				
				candidatenode = cst.child(cst.root(),seq[Occ[j].pos]);
				
				leavesprefix = cst.size(candidatenode);
					
				}
			   }
			   
			   /*find infix node*/
			   
				if(cst.is_leaf(cst.child(cst.root(),seq[Occ[j]. pos+1]))==1){
				
				candidateinfix = cst.child(cst.root(),seq[Occ[j].pos+1]);
			
				
			       }
				
			       else{
			       
			       if(cst.depth(cst.child(cst.root(),seq[Occ[j].pos+1]))<(Occ[j].size-1)){	
			       
			       auto newnode = cst.root();
				
			       auto new2node = cst.root();
			       
			       do{
								
					new2node = cst.child(newnode,(seq[Occ[j].pos+1+cst.depth(newnode)]));
					
					newnode = new2node;
					
					}
			       while(cst.depth(newnode)<(Occ[j].size-1));
			       
			       if(cst.is_leaf(newnode)==1){
				
				if(cst.depth(newnode)>(Occ[j].size-1)){
					
					candidateinfix = newnode;
					
					leavesinfix = cst.size(candidateinfix);
				
					}
				
				}
				
				else{
					
					candidateinfix = newnode;
					
					leavesinfix = cst.size(candidateinfix);
					
				}	
			       
			       }
			       else{
				
				candidateinfix = cst.child(cst.root(),seq[Occ[j].pos+1]);
				
				leavesinfix = cst.size(candidateinfix);
				
				
				}
			   }
			   
			  /*find suffix node*/  
			  
			  if(cst.is_leaf(candidateinfix)==1){
				
				candidatesuffix = candidateinfix;
				leavessuffix = leavesinfix;	
				
			  }
			  else{
				
				if(cst.depth(candidateinfix)==(Occ[j].size-1)){
				  candidatesuffix = cst.child(candidateinfix, (const char) Occ[j]. letter);
				  leavessuffix = cst.size(candidatesuffix);
				  
				}
				else{
					candidatesuffix = candidateinfix;
				  leavessuffix = leavesinfix;
				 
					}       
			 }

			 /*compute std*/
			 
				long double absentmax = 0.0;
			
				long double absentsign = sqrt((leavesprefix*leavessuffix)/leavesinfix);
			
				long double absentnumstd = 0.0;
			
				if(absentsign > 1){absentmax=absentsign;}else{absentmax=1;}	
				
				      long double absentresult = 0-((leavesprefix*leavessuffix)/leavesinfix);
				
				      if(absentresult==0){absentnumstd=0;}
					
					    else{absentnumstd = absentresult/absentmax;}
						
					    if(numt >= absentnumstd){		
					
					    memcpy( &maw[0], &seq[Occ[j]. pos], Occ[j] . size );	    
					    maw[Occ[j] . size] = (char) Occ[j]. letter;		
					    maw[Occ[j] . size + 1] = '\0';	
					    fprintf( out_fd, "%s....", maw );
					    //fprintf( out_fd, " <%lld, %lld, %c>....", Occ[j]. pos, Occ[j] .size, (char) Occ[j]. letter);

				fprintf ( out_fd, "std: %LF\n", absentnumstd );	 
				
				}           
			   
			free ( maw );      	
			}
			end = gettime();
			fprintf( stderr, " Absent Avoided Words computation: %lf secs\n", end - start);
			free ( Occ );
		}
		
	}  

        fprintf ( out_fd, ".............................................\n");  
        
	if ( fclose ( out_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}

	remove( INPUT_STR );

 	return ( 1 );       	
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
	
	
inline void compute_frequencynumk(const node_type &v, INT nuk, unsigned char * seq){
	
	auto Node = v;
	
	numk = nuk;
	
	DFSnum1 = 0;	
	
	for(auto numdegree = cst.degree(Node); numdegree>0; numdegree--){
		
		   auto NodeChild = cst.select_child(Node, numdegree);
		   
		   DFSnum1++;
		       
		   DFSunvisited1[DFSnum1]=NodeChild;
		   
	}
	
	while(DFSnum1!=0){
		
		auto DFSv1 = DFSunvisited1[DFSnum1];
			
		DFSnum1--;
		
		if(cst.is_leaf(DFSv1)==1){
			
			if(numk<(cst.depth(DFSv1))&&numk>(cst.depth(cst.parent(DFSv1)))){find_frequency(DFSv1, seq);}
			
			}
			
			else{
				
				if(numk<=(cst.depth(DFSv1))&&numk>(cst.depth(cst.parent(DFSv1)))){find_frequency(DFSv1, seq);}
					
					else{
						
						for(auto numdegree = cst.degree(DFSv1); numdegree>0; numdegree--){
		
		           auto NodeChild = cst.select_child(DFSv1, numdegree);
		
		           DFSnum1++;
		       
		           DFSunvisited1[DFSnum1]=NodeChild;}
						
						}
				}	  
}		
}
	
inline void find_frequency(const node_type &v, unsigned char * seq){	
	
	
	long double countnodenow = nodenow(v);
	
	long double countnodesuffix = nodesuffix(v);
	
	long double countnodeprefix = nodeprefix(v);
	
	long double countnodeinfix = nodeinfix(v);
	
	long double max;
	
	long double sign = sqrt((countnodeprefix*countnodesuffix)/countnodeinfix);
	
	numstd1 = 0.0;
	
	if(sign > 1){max=sign;}else{max=1;}	
		
		long double result = countnodenow-((countnodeprefix*countnodesuffix)/countnodeinfix);
		
		if(result==0){numstd=0;}
			
			else{numstd1 = result/max;}	
				
		 if(numt1 <= numstd1){	
		 	
		 	INT findposchar=cst.sn(cst.leftmost_leaf(v));
	
	     for(INT f=0; f<numk;f++){coutfd1[f]=seq[f+findposchar];}		
	      
	     fprintf ( out_fd, "%s....", (char *) coutfd1 );
	      
	     fprintf ( out_fd, "std: %LF\n", numstd1 );	    
	  
	}	
	
	
}

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

unsigned int LCParray ( unsigned char *text, INT n, INT * SA, INT * ISA, INT * LCP )
{										
	INT i=0, j=0;

	LCP[0] = 0;
	for ( i = 0; i < n; i++ ) // compute LCP[ISA[i]]
		if ( ISA[i] != 0 ) 
		{
			if ( i == 0) j = 0;
			else j = (LCP[ISA[i-1]] >= 2) ? LCP[ISA[i-1]]-1 : 0;
			while ( text[i+j] == text[SA[ISA[i]-1]+j] )
				j++;
			LCP[ISA[i]] = j;
		}

	return ( 1 );
}

unsigned int compute_maw ( unsigned char * seq, unsigned char * seq_id, struct TSwitch sw, TMaw ** Occ, INT * NOcc )
{
	INT * SA;
	INT * LCP;
	INT * invSA;
	INT n = strlen ( ( char * ) seq );
	int sigma;
    	bit_vector * Before;
    	bit_vector * Beforelcp;

	if      ( ! strcmp ( "DNA", sw . alphabet ) )   sigma = strlen ( ( char * ) DNA );
        else if ( ! strcmp ( "PROT", sw . alphabet ) )  sigma = strlen ( ( char * ) PROT );

        /* Compute the suffix array */
        SA = ( INT * ) malloc( ( n ) * sizeof( INT ) );
        if( SA == NULL )
        {
                fprintf(stderr, " Error: Cannot allocate memory for SA.\n" );
                return ( 0 );
        }

	#ifdef _USE_64
        if( divsufsort64( seq, SA,  n ) != 0 )
        {
                fprintf(stderr, " Error: SA computation failed.\n" );
                exit( EXIT_FAILURE );
        }
	#endif

	#ifdef _USE_32
        if( divsufsort( seq, SA,  n ) != 0 )
        {
                fprintf(stderr, " Error: SA computation failed.\n" );
                exit( EXIT_FAILURE );
        }
	#endif

        /*Compute the inverse SA array */
        invSA = ( INT * ) calloc( n , sizeof( INT ) );
        if( invSA == NULL )
        {
                fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
                return ( 0 );
        }

        for ( INT i = 0; i < n; i ++ )
        {
                invSA [SA[i]] = i;
        }

	LCP = ( INT * ) calloc  ( n, sizeof( INT ) );
        if( LCP == NULL )
        {
                fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
                return ( 0 );
        }

        /* Compute the LCP array */
        if( LCParray( seq, n, SA, invSA, LCP ) != 1 )
        {
                fprintf(stderr, " Error: LCP computation failed.\n" );
                exit( EXIT_FAILURE );
        }

	free ( invSA );

	INT v_size = 2 * n;

    	Before = new bit_vector[sigma];
        if( Before == NULL )
        {
                fprintf(stderr, " Error: Cannot allocate memory for Before.\n" );
                return ( 0 );
        }
	Beforelcp = new bit_vector[sigma];
        if( Beforelcp == NULL )
        {
                fprintf(stderr, " Error: Cannot allocate memory for BeforeLCP.\n" );
                return ( 0 );
        }
    	for ( INT i = 0; i < sigma; i++ )
	{
		Before[i] = bit_vector( v_size, 0 );
		Beforelcp[i] = bit_vector( v_size, 0 );
    	}

	GetBefore ( seq, n , sigma, SA, LCP, Before, Beforelcp );

	GetMaws( seq, seq_id, SA, n, sigma, LCP, Before, Beforelcp, sw . k, sw . K, sw . output_filename, Occ, NOcc );

	delete [] Before;
	delete [] Beforelcp ;
	free ( SA );
	free ( LCP );

	return ( 1 );
}


unsigned char Mapping( int a )
{
	char c = DEL;
        switch ( a )
	{
            case 0:
                c = 'A';
                break;
            case 1:
                c = 'C';
                break;
            case 2:
                c = 'G';
                break;
            case 3:
                c = 'T';
                break;
            case 4:
                c = 'N';
                break;
            case 5:
                c = 'R';
                break;
            case 6:
                c = 'D';
                break;
            case 7:
                c = 'Q';
                break;
            case 8:
                c = 'E';
                break;
            case 9:
                c = 'H';
                break;
            case 10:
                c = 'I';
                break;
            case 11:
                c = 'L';
                break;
            case 12:
                c = 'K';
                break;
            case 13:
                c = 'M';
                break;
            case 14:
                c = 'F';
                break;
            case 15:
                c = 'P';
                break;
            case 16:
                c = 'S';
                break;
            case 17:
                c = 'W';
                break;
            case 18:
                c = 'Y';
                break;
            case 19:
                c = 'V';
                break;
        }
	return ( c );
}

int RevMapping ( unsigned char b )
{
	int a = -1;
        switch ( b )
	{
            case 'A':
                a = 0;
                break;
            case 'C':
                a = 1;
                break;
            case 'G':
                a = 2;
                break;
            case 'T':
                a = 3;
                break;
            case 'N':
                a = 4;
                break;
            case 'R':
                a = 5;
                break;
            case 'D':
                a = 6;
                break;
            case 'Q':
                a = 7;
                break;
            case 'E':
                a = 8;
                break;
            case 'H':
                a = 9;
                break;
            case 'I':
                a = 10;
                break;
            case 'L':
                a = 11;
                break;
            case 'K':
                a = 12;
                break;
            case 'M':
                a = 13;
                break;
            case 'F':
                a = 14;
                break;
            case 'P':
                a = 15;
                break;
            case 'S':
                a = 16;
                break;
            case 'W':
                a = 17;
                break;
            case 'Y':
                a = 18;
                break;
            case 'V':
                a = 19;
                break;
        }
	return ( a );
}

unsigned int GetBefore (
				unsigned char * seq,
				INT n,
				int sigma,
				INT * SA,
				INT * LCP,
				bit_vector * Before,
				bit_vector * Beforelcp )
{
        INT hm = 0;
        INT k = 0;
        INT lcp;
        INT mem;
        INT proxa;
        INT proxb;

        TStack lifo_lcp;
        StackNew ( &lifo_lcp, sizeof( INT ) );
        TStack lifo_mem;
        StackNew ( &lifo_mem, sizeof( INT ) );
        TStack lifo_rem;
        StackNew ( &lifo_rem, sizeof( INT ) );

        lcp = 0;
        StackPush(&lifo_lcp, &lcp);

        /* Max LCP value */
        for ( int i = 0; i < n; i++ )
                if( LCP[i] > hm )
                        hm = LCP[i];
        hm = hm + 2;

        bit_vector* interval = new bit_vector[sigma];
        for ( INT i = 0; i < sigma; i++)
                interval[i]=bit_vector(hm,0);

	      interval[RevMapping(seq[n- 1])][0]=1;

        // First pass : top-down
        for ( INT i = 0; i < n; i++ )
        {
                // first we update the interval table
                // we empty the interval that corresponds to a higher lcp value
                if ( i > 0 && LCP[i] < LCP[i-1])
                {
                        StackPop(&lifo_lcp,&lcp);
                        while(!StackEmpty(&lifo_lcp)&&lcp>LCP[i])
                        {
                                StackPop(&lifo_lcp,&mem);
                                if (mem <=LCP[i])
                                {
                                        for (int j=0; j<sigma; j++)
                                        {
                                                if (mem!=LCP[i]){interval[j][LCP[i]]=interval[j][lcp];}     //initialisation of the next intervals if it hasn't been open
                                                Before[j][2*i-1]=interval[j][lcp];
                                                Beforelcp[j][2*i-1]=interval[j][lcp];
                                                if (mem==LCP[i]){Beforelcp[j][2*i-1]=interval[j][mem];}
                                        }
                                }

                                for (int j =0; j<sigma ;j++){interval[j][lcp]=0;}
                                        lcp=mem;
                        }
                        StackPush(&lifo_lcp,&lcp);
                }

        	// we update those having a lower lcp
                if ( SA[i] - 1 >= 0 )
                        k = RevMapping(seq[SA[i] - 1]);
                else
                        k = - 1;
                if ( k != -1 )
                {
                        while(!StackEmpty(&lifo_lcp))
                        {
                                StackPop(&lifo_lcp,&lcp);
                                StackPush(&lifo_mem, &lcp);
                                if (interval[k][lcp]==1){break;}
                                interval[k][lcp]=1;
                        }

                        while (!StackEmpty(&lifo_mem))
                        {
                                StackPop(&lifo_mem,&lcp);
                                StackPush(&lifo_lcp, &lcp);
                        }

                        interval[k][LCP[i]]=1;
                }

                if ( ( i - 1 ) >= 0 && SA[i - 1] - 1 >= 0 )
                        if (i>0 && LCP[i]>0 && RevMapping(seq[SA[i - 1] - 1])!=-1) // in this case we also add the letter preceding the last suffix
                        {
                                interval[RevMapping(seq[SA[i - 1] - 1])][LCP[i]]=1;
                        }

                // Before, Before_LCP
                for (int j =0; j<sigma ;j++)
                {
                        Beforelcp[j][2*i]=interval[j][LCP[i]];
                }

                if (k!=-1)
                {
                        Before[k][2*i+1]=1;
                        Before[k][2*i]=1;
                        Beforelcp[k][2*i+1]=1;
                        Beforelcp[k][2*i]=1;
                }

                StackPop(&lifo_lcp,&lcp);
                if(lcp!=LCP[i])  // no duplicates
                {
                        StackPush(&lifo_lcp,&lcp);
                        lcp=LCP[i];
                }
                StackPush(&lifo_lcp,&lcp);
        }

        //second pass : bottom-up
        //we empty the interval table

        while(!StackEmpty(&lifo_lcp))
        {
                StackPop(&lifo_lcp,&lcp);
                for (int j =0; j<sigma ;j++){interval[j][lcp]=0;}
        }
        lcp=0;
        StackPush(&lifo_lcp,&lcp);

        for (INT i=n-1; i>-1; i--)
        {
                StackPop(&lifo_lcp,&lcp);
                proxa=LCP[i]+1;   //proxb is the lcp-value that is just higher than LCP[i]
                while(!StackEmpty(&lifo_lcp) && lcp>LCP[i])
                {
                        StackPush(&lifo_rem,&lcp);
                        StackPop(&lifo_lcp,&mem);
                        if (mem <LCP[i])            //initialisation of the interval if it hasn't been open
                        {
                                for (int j=0; j<sigma; j++){interval[j][LCP[i]]=interval[j][lcp];}
                                proxa=lcp;
                        }
                        if (mem ==LCP[i]){proxa=lcp;}
                        lcp=mem;
                }
                StackPush(&lifo_lcp,&lcp);

                // we update the lower intervals
                for (int k=0; k<sigma;k++)
                {
                        if(Before[k][2*i]==1)
                        {
                                while(!StackEmpty(&lifo_lcp))
                                {
                                        StackPop(&lifo_lcp,&lcp);
                                        StackPush(&lifo_mem, &lcp);
                                        if (interval[k][lcp]==1){break;}
                                        interval[k][lcp]=1;
                                }

                                while (!StackEmpty(&lifo_mem))
                                {
                                        StackPop(&lifo_mem,&lcp);
                                        StackPush(&lifo_lcp, &lcp);
                                }
                                interval[k][LCP[i]]=1;
                        }
                }

                for (int j =0; j<sigma ;j++)
                { //update interval and Before
                        Beforelcp[j][2*i]=Beforelcp[j][2*i] || interval[j][LCP[i]];

                        if(i<n-1)
                        {
                                Before[j][2*i+1]=Before[j][2*i+1] || interval[j][proxb];  //proxb is the lcp-value that is just higher than LCP[i+1]
                                Beforelcp[j][2*i+1]=interval[j][LCP[i+1]] || Beforelcp[j][2*i+1];
                        }
                }

                proxb=proxa;

                //we suppress higher intervals
                if (i<n-1 && LCP[i+1]>LCP[i])
                {
                        StackPop(&lifo_rem,&lcp);  // this lcp is the one that is just higher than LCP[i]
                        for (int j =0; j<sigma ;j++)
                        {
                                Before[j][2*i]=Before[j][2*i] || interval[j][lcp];
                                interval[j][lcp]=0;
                        }

                        while(!StackEmpty(&lifo_rem))
                        {
                                StackPop(&lifo_rem,&lcp);
                                for (int j =0; j<sigma ;j++){interval[j][lcp]=0;}

                        }
                }

                StackPop(&lifo_lcp,&lcp);
                if(lcp!=LCP[i])  // no duplicates
                {
                        StackPush(&lifo_lcp,&lcp);
                        lcp=LCP[i];
                }
                StackPush(&lifo_lcp,&lcp);
        }

        delete[] interval;
        StackDispose(&lifo_lcp);
        StackDispose(&lifo_mem);
        StackDispose(&lifo_rem);

	return ( 1 );
}

unsigned int GetMaws( unsigned char * seq, unsigned char * seq_id, INT * SA, INT n, int sigma, INT * LCP, bit_vector* Before, bit_vector* Beforelcp, unsigned int k, unsigned int K, char * out_file, TMaw ** Occ, INT * NOcc )
{
    	FILE * out_fd;
	char * maw;

	// compute a bitvector that contains a `1', if an identical row has already been seen => to avoid duplicates.
    	bit_vector Mem = bit_vector(n,0);
    	TStack lifo_lcp;
	StackNew (&lifo_lcp, sizeof( INT ) );
	INT lcp = 0;
	INT mem;
	StackPush(&lifo_lcp, &lcp);

	for ( INT i = 0; i < n; i++ )
	{
        	StackPop(&lifo_lcp,&lcp);
            	while(!StackEmpty(&lifo_lcp)&&lcp>LCP[i])
            	{
                	StackPop(&lifo_lcp,&mem);
                	if ( mem == LCP[i] )
                	{
                    		Mem[i] = 1;
                	}
                	lcp = mem;
            	}
            	StackPush(&lifo_lcp,&lcp);
            	lcp = LCP[i];
            	StackPush(&lifo_lcp,&lcp);
	}
	StackDispose(&lifo_lcp);

        maw = ( char * ) calloc( ( K + 1 ) , sizeof( char ) );
        if( maw == NULL )
        {
                fprintf(stderr, " Error: Cannot allocate memory.\n" );
                return ( 0 );
        }
        
	vector<mytuple> xtuple;
            
        for ( INT i = 0; i < n; i++ )
    	{
        	for( int l = 0; l < sigma; l++ )
        	{
            		bool B1 = (
					Before[l][2 * i] == 0 &&
					Beforelcp[l][2 * i] == 1 &&
					SA[i] + LCP[i] < n &&
					LCP[i] + 2 <= K	&&
					LCP[i] + 2 >= k	);

            		bool B2 = (
					i < n - 1 &&
					Before[l][2 * i + 1] == 0 &&
					Beforelcp[l][2 * i + 1] == 1 &&
					SA[i] + LCP[i + 1] < n &&
					LCP[i + 1] + 2 <= K &&
					LCP[i + 1] + 2 >= k &&
					Mem[i + 1] == 0 );
			
			/* Here we report MAWs */
		    	if ( B1 )
		    	{
				maw[0] = Mapping( l );
				INT start = SA[i];
				INT size = SA[i]+ LCP[i] + 1 - start;
				memcpy( &maw[1], &seq[start], size );
				maw[size + 1] = '\0';
				xtuple . push_back(make_tuple(( INT ) maw[0], i, start, size ));

		    	}
		    	else if ( B2 )
		    	{
				maw[0] = Mapping( l );
				INT start = SA[i];
				INT size = SA[i] + LCP[i + 1] + 1 - start;
				memcpy( &maw[1], &seq[start], size );
				maw[size + 1] = '\0';
				xtuple . push_back(make_tuple(( INT ) maw[0], i, start, size ));
		    	}
        	}

    	}
    	
    	std::sort( xtuple.begin(), xtuple.end() );
	
	for(vector<mytuple>::iterator iter = xtuple.begin(); iter != xtuple.end(); iter++)
	  ( * NOcc ) ++;
	
	( * Occ ) = ( TMaw * ) realloc ( ( * Occ ), ( ( * NOcc ) ) * sizeof ( TMaw ) );
	
	INT j = 0;
	for(vector<mytuple>::iterator iter = xtuple.begin(); iter != xtuple.end(); iter++)
	{
		( *Occ )[j] . letter = get<0>(*iter);
		( *Occ )[j] . pos = get<2>(*iter);
		( *Occ )[j] . size = get<3>(*iter);
		j++;
  	}

	free ( maw );

	return ( 1 );
}

////////////////////////////////////////////////////////////////////

inline void compute_all_of_occuring_avoided_words(const node_type &v, unsigned char * seq){
	
	auto Node = v;	
	auto SuffixNode=cst.root();
	
	INT DFSall = 0;
	INT findposchar1;		
	INT f1;	
	INT countleafprefix=0;
  INT countleafcurrent=0;
  INT countleafinfix=0;
  INT countleafsuffix=0;
  
  long double countnodenow1;
  long double countnodesuffix1;
	long double countnodeprefix1;
	long double countnodeinfix1; 	
	long double max1 = 0.0;
  long double sign1 = 0.0;
  long double result1 = 0.0;
  long double numstd1 = 0.0;
	
	for(auto num_degree_all_1 = cst.degree(Node); num_degree_all_1>0; num_degree_all_1--){
		
	    auto NodeChild1 = cst.select_child(Node, num_degree_all_1);
	    
	    if((cst.is_leaf(NodeChild1)==0)&&(cst.depth(NodeChild1)>1)){
			
			    DFSall++;
		       
		      DFSunvisitedall[DFSall]=NodeChild1; 
	                
	    }
	    
	    if((cst.is_leaf(NodeChild1)==0)&&(cst.depth(NodeChild1)==1)){
	    	
	    	  for(auto num_degree_all_2 = cst.degree(NodeChild1); num_degree_all_2>0; num_degree_all_2--){
		
		          auto NodeChild2 = cst.select_child(NodeChild1, num_degree_all_2);
		          
		          if(cst.is_leaf(NodeChild2)==0){
			
			            DFSall++;
		       
		              DFSunvisitedall[DFSall]=NodeChild2; 
		              
		          }
	                
	        }
	            
	    }
	         
  }
	   
	            
	while(DFSall!=0){
		
	     auto DFSallnode = DFSunvisitedall[DFSall];
			
	     DFSall--;
			
			 auto PrefixNode = DFSallnode;
			 
			 countleafprefix = 0;
			 
		   countleafprefix = cst.size(PrefixNode);
		   
		   auto InfixNode = cst.sl(PrefixNode);
		   
		   countleafinfix = 0;
		   
		   if(InfixNode==cst.root()){countleafinfix = 0;}
		   	
		   	else{countleafinfix = cst.size(InfixNode);}
			
			 for(auto num_degree_all_3 = cst.degree(PrefixNode); num_degree_all_3>0; num_degree_all_3--){
		
		       auto NodeChild3 = cst.select_child(PrefixNode, num_degree_all_3);
		   
		       auto CurrentNode = NodeChild3;
		        
		       if((cst.is_leaf(CurrentNode)==1)&&(((cst.depth(CurrentNode))-(cst.depth(PrefixNode)))==1)){countleafcurrent = 0;}
		              	
		       else{
		        
		           countleafcurrent = 0;
		           
		           if(cst.is_leaf(CurrentNode)==1){countleafcurrent=1;}
		           	else{countleafcurrent = cst.size(CurrentNode);}
		           		
	//method 1 to find suffix node:
		           		
		           		auto CandidateSuffixNode = cst.sl(CurrentNode);
		           		
		           		while(cst.parent(CandidateSuffixNode)!=InfixNode){CandidateSuffixNode=cst.parent(CandidateSuffixNode);}
		           		
		           		SuffixNode = CandidateSuffixNode;

		           
  //method 2 to find suffix node:
  
              // auto CharCurrent = cst.edge(CurrentNode, ((cst.depth(PrefixNode))+1));          
		           
		          // for(auto num_degree_all_4 = cst.degree(InfixNode); num_degree_all_4>0; num_degree_all_4--){
		
		          //    auto NodeChild4 = cst.select_child(InfixNode, num_degree_all_4);
		  
		          //    auto CandidateSuffixNode = NodeChild4;
		              
		          //    if(cst.edge(CandidateSuffixNode, ((cst.depth(InfixNode))+1))==CharCurrent){SuffixNode=CandidateSuffixNode;}
		           
		          // }
		           
  //method 3 to find suffi xnode: 
 
                // auto CharCurrent = cst.edge(CurrentNode, ((cst.depth(PrefixNode))+1));
		  
		            // SuffixNode = cst.child(InfixNode, CharCurrent);
		            
		            //
		        
		           countleafsuffix = 0;
		           
		           if(SuffixNode==cst.root()){countleafsuffix = 0;}
		           
		           else if(cst.is_leaf(SuffixNode)==1){countleafsuffix=1;}
		           	
		           	else if(cst.is_leaf(SuffixNode)==0){countleafsuffix = cst.size(SuffixNode);}
		           
		           /////////formula//////////

		                  countnodenow1 = countleafcurrent;
	
                      countnodesuffix1 = countleafsuffix;
	
	                    countnodeprefix1 = countleafprefix;
	
	                    countnodeinfix1 = countleafinfix;
	
	                    max1 = 0.0;
	
          	          sign1 = 0.0;
          	          
          	          sign1 = (long double)sqrt((countnodeprefix1*countnodesuffix1)/countnodeinfix1);
	
          	          numstd1 = 0.0;
	
          	          if(sign1 > 1.0){max1=sign1;}else{max1=1.0;}	
		
	        	          result1 = 0.0;
	        	          
	        	          result1 = (long double)(countnodenow1-(long double)((countnodeprefix1*countnodesuffix1)/countnodeinfix1));
		
	        	          if(result1==0.0){numstd1=0.0;}
			
	      		          else{numstd1 = ((long double)(result1/max1));}	
				
		                  if(numt >= numstd1){	
		                  	
		                  	  findposchar1=0;
		 	
		      	              findposchar1=cst.sn(cst.leftmost_leaf(CurrentNode));

	                        for(f1=0; f1<((cst.depth(PrefixNode))+1);f1++){coutfdall[f1]=seq[f1+findposchar1];}		
	      
	                        fprintf ( out_fd, "%s....", (char *) coutfdall);
	      
	                        fprintf ( out_fd, "std: %LF\n", numstd1 );
	                        
	                        memset((char *)coutfdall, 0, nnnn );
 
	  
	                    }	
	                    
	             //
	                    
	             if(cst.is_leaf(CurrentNode)==0){
		
		               DFSall++;
		       
		               DFSunvisitedall[DFSall]=CurrentNode;	
		                  
		           }
		              	
		        }
		              
		    } 
		       
		}
		  
}
		              
