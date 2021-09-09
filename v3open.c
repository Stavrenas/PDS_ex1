#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <sys/time.h>
#include "mmio.h"
#include <omp.h>



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
        col_l = col_coo[l] - isOneBased;

        uint32_t dst = col[col_l];
        row[dst] = row_coo[l]+1;

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


    int exists( uint32_t * const row, uint32_t * const col,uint32_t ifind, uint32_t jfind)//function to check if the element is nonzero
    {
        uint32_t col_elements=col[jfind]-col[jfind-1];
		int ret=0;
		int isOneBased=1; 
        if(col_elements>0)
        {
            for(uint32_t i=0; i <col_elements; i++)
            {
				//printf("Checking [%d][%d]"
                if(row[col[jfind-1]+i]==ifind)
                    ret=1;
            }

        }
		
		return ret;

    }
	
	
	void V3open( uint32_t * const row, uint32_t * const col, uint32_t  n){
		
		uint32_t *triangles= (uint32_t*)calloc(n, sizeof(uint32_t)) ;

		#pragma omp parallel 
		{
		//printf("thread %d, of %d\n",omp_get_thread_num(),omp_get_num_threads());
		int nthreads=omp_get_num_threads();
		int id=omp_get_thread_num();

        for(uint32_t i=id+1; i<=n; i+=nthreads){
			
			

			for(uint32_t i1=col[i-1]; i1<col[i]; i1++) // iterate through elements in jth row
			{
                   if(row[i1]==i)
					   continue;
				   uint32_t j=row[i1];
				   uint32_t start=col[j-1];
				   
			       for (uint32_t i2=col[i-1]; i2<col[i] ; i2++){
					   
						for (uint32_t i3=start; i3<col[j] && row[i3]<=row[i2];i3++){

							if(row[i3]==row[i2]){

							
							#pragma omp critical
							{
							triangles[j-1]++;
							triangles[row[i2]-1]++;
							triangles[i-1]++;
							start=i3;
							}
							

							
							
							break;                               //do not check the rest


					}
					
				} 
			}
     	}
		
		}
		}
		
    uint32_t sum=0;
		
    for(uint32_t i=0; i<n; i++)
    {
        //printf( "i: %d ->%ld \n",i+1,triangles[i]);
        sum+=triangles[i];
    }

    printf("TrianglesV3open: %d\n",sum/3);
	free(triangles);
		
	}
	
	
	
	void V3open2( uint32_t * const row, uint32_t * const col, uint32_t  n){
		
		omp_set_num_threads(4); 
		int nthreads=omp_get_num_threads();
		uint32_t **triangles= (uint32_t**)calloc(n, sizeof(uint32_t*)) ;
		for(int j=0; j<n; j++){
			triangles[j]=(uint32_t*)calloc(nthreads, sizeof(uint32_t));
		 }
		
		int i;
		#pragma omp parallel shared (i)
		{
		int id=omp_get_thread_num()+1;

        for(uint32_t i=1; i<=n-2; i++){
			if(i%id==0)
			{
			
			for(uint32_t i1=col[i-1]; i1<col[i]; i1++) // iterate through elements in jth row
			{
                   if(row[i1]==i)
					   continue;
				   uint32_t j=row[i1];
				   uint32_t start=col[j-1];
				   
			       for (uint32_t i2=col[i-1]; i2<col[i] ; i2++){
					   
						for (uint32_t i3=start; i3<col[j] && row[i3]<=row[i2];i3++){

							if(row[i3]==row[i2]){

						
							
							triangles[j-1][id-1]++;
							triangles[row[i2]-1][id-1]++;
							triangles[i-1][id-1]++;
							start=i3;
							
							

							
							
							break;                               //do not check the rest


					}
					
				} 
			}
     	}
		
		}
		}
	}
		
    uint32_t sum=0;
		
    for(uint32_t i=0; i<n; i++)
    {
        //printf( "i: %d ->%ld \n",i+1,triangles[i]);
		for(int j=0;j<nthreads;j++)
          sum+=triangles[i][j];
    }

    printf("TrianglesV3open: %d\n",sum/3);
	free(triangles);
		
	}
	



    int main(int argc, char *argv[])
    {

      int ret_code;
    MM_typecode matcode;
    FILE *f;
    uint32_t M, N, nz;   
    uint32_t i, *I, *J;
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

    I = (uint32_t *) malloc(nz * sizeof(uint32_t));
    J = (uint32_t *) malloc(nz * sizeof(uint32_t));
    val = (double *) malloc(nz * sizeof(double));


    /* NOTE: when reading in doubles, ANSI C requires the use of the "l"  */
    /*   specifier as in "%lg", "%lf", "%le", otherwise errors will occur */
    /*  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            */

    for (i=0; i<nz; i++)
    {
        fscanf(f, "%d %d", &I[i], &J[i]);
        I[i]--;  /* adjust from 1-based to 0-based */
        J[i]--;
    }

    if (f !=stdin) fclose(f);

    /************************/
    /* now write out matrix */
    /************************/


        //mm_write_banner(stdout, matcode);
        //mm_write_mtx_crd_size(stdout, M, N, nz);
        //Up to this point, I[] cointains row index and J[] column index for the nonzero elements
        /*for (i=0; i<nz; i++)
        {
             fprintf(stdout, "%d %d\n", I[i]+1, J[i]+1);
             printf("%d %d\n", I[i]+1, J[i]+1);
        }
		*/


        const uint32_t nnz = nz;
        const uint32_t n   = M;  //from  mm_read_mtx_crd_size(f, &M, &N, &nz)
		printf("M is %d, nnz is %d\n",n,nnz);
        uint32_t * csc_row = (uint32_t *)malloc(nnz     * sizeof(uint32_t));
        uint32_t * csc_col = (uint32_t *)malloc((n + 1) * sizeof(uint32_t));
        uint32_t isOneBased = 0;
        // Call coo2csc for isOneBase false
        coo2csc(csc_row, csc_col, I, J, nnz, n, isOneBased);

        struct timeval tStart ;
        for(int i=0; i<5;i++){
		tStart = tic();
		
        V3open(csc_row,csc_col,n);
		
		printf("%.6f\n", toc(tStart));
		}

        /* cleanup variables */
        free( csc_row );
        free( csc_col );

    }
