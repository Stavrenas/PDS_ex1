#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <sys/time.h>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <pthread.h> //pthread library
#include "mmio.h"



struct timeval tic()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv;
}

double toc(struct timeval begin)
{
    struct timeval end;
    gettimeofday(&end, NULL);
    double stime = ((double) (end.tv_sec - begin.tv_sec) * 1000 ) +
                   ((double) (end.tv_usec - begin.tv_usec) / 1000 );
    stime = stime / 1000;
    return(stime);
}




void coo2csc(
    uint32_t       * const row,       /*!< CSC row start indices */
    uint32_t       * const col,       /*!< CSC column indices */
    uint32_t const * const row_coo,   /*!< COO row indices */
    uint32_t const * const col_coo,   /*!< COO column indices */
    uint32_t const         nnz,       /*!< Number of nonzero elements */
    uint32_t const         n,         /*!< Number of rows/columns */
    uint32_t const         isOneBased /*!< Whether COO is 0- or 1-based */
)
{

    for (uint32_t l = 0; l < n+1; l++) col[l] = 0;


    for (uint32_t l = 0; l < nnz; l++)
        col[col_coo[l] - isOneBased]++;

    // ----- cumulative sum
    for (uint32_t i = 0, cumsum = 0; i < n; i++)
    {
        uint32_t temp = col[i];
        col[i] = cumsum;
        cumsum += temp;
    }
    col[n] = nnz;
    // ----- copy the row indices to the correct place
    for (uint32_t l = 0; l < nnz; l++)
    {
        uint32_t col_l;
	    /*if(col_coo[l]==row_coo[l])
	    continue;
	    */
        col_l = col_coo[l] - isOneBased;

        uint32_t dst = col[col_l];
        row[dst] = row_coo[l];   //row[dst] = row_coo[l]-isOneBased; this caused segmentation fault

        col[col_l]++;
    }
    // ----- revert the column pointers
    for (uint32_t i = 0, last = 0; i < n; i++)
    {
        uint32_t temp = col[i];
        col[i] = last;
        last = temp;
    }

}




	
	
	int binsearch (uint32_t * const row, uint32_t * const col,uint32_t ifind, uint32_t start, uint32_t elements){
		uint32_t middle=elements/2 + start; 
		if(middle<start)
			return 0;
		if(row[middle]==ifind){
		      return 1;
		}
		else if((elements<=1 && row[middle]!=ifind )){
             return 0;
		}		
		else if(row[middle] < ifind ){
			binsearch(row, col, ifind, middle, (elements - elements/2));
		}
		else{
			binsearch(row, col, ifind, start , (elements/2));
		}
		return -1;
	}
		
	
	int search( uint32_t * const row, uint32_t * const col,uint32_t ifind, uint32_t jfind)//function to check if the element is nonzero
    {
		if(ifind==jfind){
			return 0;  //since we do not care for self edges
		}
		if(ifind<jfind){
			uint32_t temp=ifind;
			ifind=jfind;
			jfind=ifind;
		}	
        uint32_t elements=col[jfind]-col[jfind-1];
		int start=col[jfind-1];
		int middle;
		//binsearch(row,col,ifind,start,elements);
		middle=start+elements/2;
		while (elements >0){
			if(row[middle]==ifind)
                 return 1;
            if(row[middle]!=ifind && elements==1)
                 return 0;				 
			if(row[middle]>ifind){
				elements-=elements/2;
				if (elements%2==1)
				  middle-=elements/2-1;
			    else
				  middle-=elements/2;
			}
			else{
				if(elements%2==0)
					elements-=elements/2-1;
				else
			        elements-=elements/2;
				middle+=elements/2;		
			}
			
		}
		return -1;

    }
	

	
	
	/*If I multiply the sparce binary summetric matrix A with the dense vector e that
	has each element set to 1, then I get
	C=Ae 
	C[1]=A11*e1+A12*e2+A13*e3+...= sum of A1j for j:1-n since ei=1
	This means that the vector C is equal to the number of nnz elements of 
	row i in matrix A. This is equal to the vector that contains the number of 
	elements in each row in csc format, beacause Aij=0 or 1.

    If A is not binary, then C[i]=sum of the values of the nonzero elements in i_th row     
	
	
	
	
	let A*A=B. B[i][j] = ith_row * jth_column == ith_row*jth_row which equals to the 
	sum of elements of i_th and j_th row that belong in the same column

	*/
	
	
	
int find_common( uint32_t * const row, uint32_t * const col,uint32_t i, uint32_t j){
	int start=col[j-1];
	uint32_t t=0;
	//uint32_t *t= (uint32_t*)calloc(__cilkrts_get_nworkers(), sizeof(uint32_t)) ;
	// __cilkrts_get_worker_number() and __cilkrts_get_nworkers()
	
	
				 for (uint32_t i2=col[i-1]; i2<col[i] ; i2++){  

                    for(uint32_t i3=start; i3<col[j] && row[i3]<=row[i2] ;i3++){
						
						if(row[i3]==row[i2]){
							t++;
							start=i3;
							break;//do not check the rest

						}
					} 
				 }

	return t;

}	
	
	
	
	
	void V4cilk( uint32_t * const row, uint32_t * const col, uint32_t  n){
		uint32_t *triangles= (uint32_t*)calloc(n, sizeof(uint32_t)) ;
		uint32_t *C= (uint32_t*)calloc(col[n], sizeof(uint32_t)) ; //the value array equivalent in csc format
		uint32_t j,t;


       cilk_for (uint32_t i=1; i<=n; i++){
		   

			   cilk_for (uint32_t i1=col[i-1]; i1<col[i]; i1++){	 
				  j=row[i1];
				 //C[k]=cilk_spawn find_common(row,col,i,j);
				 C[i1]= find_common(row,col,i,j);

				 
			   }	   
		}
	    
    //cilk_sync;
    uint32_t sum=0;
       for(uint32_t i=1; i<=n; i++){
			   for(uint32_t i1=col[i-1]; i1<col[i]; i1++){
				   triangles[i-1]+=C[i1]; //normally it would be C[k]/2 but in the case that C[k] is odd, then 1 element is lost
			   }
	   }
	   
	for(uint32_t i=0; i<n; i++)
    {    
        triangles[i]/=2;
        //printf("%d - %d\n",i+1,triangles[i]/2);
        sum+=triangles[i];
    }

    printf("\nTrianglesV4cilk: %d\n",sum/(3));

	free(C);
	free(triangles);
	}
	





    int main(int argc, char *argv[])
    {

    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;   
    uint32_t i, *I, *J,*Itotal,*Jtotal;
    double *val;

    if (argc < 2)
	{
		fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
		exit(1);
	}
    else    
    { 
        if ((f = fopen(argv[1], "r")) == NULL) 
            exit(1);
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }


    /*  This is how one can screen matrix types if their application */
    /*  only supports a subset of the Matrix Market data types.      */

    if (mm_is_complex(matcode) && mm_is_matrix(matcode) && 
            mm_is_sparse(matcode) )
    {
        printf("Sorry, this application does not support ");
        printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
        exit(1);
    }

    /* find out size of sparse matrix .... */

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) !=0)
        exit(1);


    /* reseve memory for matrices */

	Itotal = (uint32_t *) malloc(nz*2 * sizeof(uint32_t));
    Jtotal = (uint32_t *) malloc(nz*2 * sizeof(uint32_t));
	int val1,val2;


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    for (i=0; i<nz; i++)
    {
        fscanf(f, "%d %d", &val1, &val2);
	   
		Itotal[2*i]=val1;
		Jtotal[2*i]=val2;
		Itotal[2*i+1]=val2;
		Jtotal[2*i+1]=val1;
		
		
    }

    if (f !=stdin) fclose(f);

    /************************/
    /* now write out matrix */
    /************************/



        
        uint32_t nnz = nz*2;
        uint32_t n   = M;  //from  mm_read_mtx_crd_size(f, &M, &N, &nz)
		printf("M is %d, nnz is %d\n",n,nnz);
        uint32_t isOneBased = 1;
        uint32_t * csc_row = (uint32_t *)malloc(nnz     * sizeof(uint32_t));
        uint32_t * csc_col = (uint32_t *)malloc((n + 1) * sizeof(uint32_t));
        coo2csc(csc_row, csc_col, Itotal, Jtotal, nnz, n, isOneBased);
		

        struct timeval tStart;
		for(int i=0; i<5;i++){
		tStart = tic();
		
        V4cilk(csc_row,csc_col,n);
		
		printf("%.6f\n", toc(tStart));
		}

        /* cleanup variables */
		
		free( csc_row );
        free( csc_col);
		free( Itotal );
        free( Jtotal );


    }
