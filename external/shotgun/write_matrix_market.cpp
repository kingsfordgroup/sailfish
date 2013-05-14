/* 
*   Matrix Market I/O example program
*
*   Create a small sparse, coordinate matrix and print it out
*   in Matrix Market (v. 2.0) format to standard output.
*
*   (See http://math.nist.gov/MatrixMarket for details.)
*
*/

#include <stdio.h>
#include <stdlib.h>
#include "assert.h"
#include "mmio.h"


int write_to_file(const char * filename, int * I, int * J, double * val, int M, int N, int nz)
{
    MM_typecode matcode;                        
    int i;

    mm_initialize_typecode(&matcode);
    mm_set_matrix(&matcode);
    mm_set_coordinate(&matcode);
    mm_set_real(&matcode);
    FILE * fp = fopen(filename,"w");
    assert(fp);
    mm_write_banner(fp, matcode); 
    mm_write_mtx_crd_size(fp, M, N, nz);

    /* NOTE: matrix market files use 1-based indices, i.e. first element
      of a vector has index 1, not 0.  */

    for (i=0; i<nz; i++)
        fprintf(fp, "%d %d %10.3g\n", I[i]+1, J[i]+1, val[i]);

    fclose(fp);
    return 0;
}

