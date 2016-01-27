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
#include "stack.h"        					  // include header for stack structure

using namespace sdsl;
using namespace std;
typedef cst_sct3<> cst_t;
typedef cst_sada<> csts_t;
typedef bp_interval<> node_type;


cst_t cst;
csts_t csts;
INT nnnn;

INT DFSnum2;
unsigned int numk = 0;
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

        fprintf ( out_fd, "k = %d \n", numk );

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
        fprintf( stderr, " Occuring Avoided Words computation: %lf secs\n", end - start);

        start = gettime();
	TMaw * Occ = NULL; 
 	unsigned int NOcc = 0;
	compute_maw ( seq, seq_id, sw, &Occ, &NOcc );
        end = gettime();
        fprintf( stderr, " Absent Avoided Words computation: %lf secs\n", end - start);

	if ( fclose ( out_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}

	remove( INPUT_STR );
       	
	free(coutfd);
       
       	free(DFSunvisited2);   

	free ( Occ );

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

unsigned int compute_maw ( unsigned char * seq, unsigned char * seq_id, struct TSwitch sw, TMaw ** Occ, unsigned int * NOcc )
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
        if( ( SA == NULL) )
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
        if( ( invSA == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
                return ( 0 );
        }

        for ( INT i = 0; i < n; i ++ )
        {
                invSA [SA[i]] = i;
        }

	LCP = ( INT * ) calloc  ( n, sizeof( INT ) );
        if( ( LCP == NULL) )
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
        if( ( Before == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for Before.\n" );
                return ( 0 );
        }
	Beforelcp = new bit_vector[sigma];
        if( ( Beforelcp == NULL) )
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

unsigned int GetMaws( unsigned char * seq, unsigned char * seq_id, INT * SA, INT n, int sigma, INT * LCP, bit_vector* Before, bit_vector* Beforelcp, unsigned int k, unsigned int K, char * out_file, TMaw ** Occ, unsigned int * NOcc )
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
        if( ( maw == NULL) )
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

