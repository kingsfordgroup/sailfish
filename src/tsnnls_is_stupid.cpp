#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Sparse>
#include "tsnnls/tsnnls.h"

using std::string;
using std::vector;
using std::ifstream;

template <typename MatT>
bool readMatrix( const string& fname, MatT*& M ) {
	typedef Eigen::Triplet<double> T;
	ifstream ifile(fname);
	std::string junk;
	size_t nrows, ncols, nnz;
	ifile >> junk;
	ifile >> nrows >> ncols;
	ifile >> nnz;

	M = new MatT(nrows, ncols);

	std::vector<T> triplets;
	triplets.reserve(nnz);
	size_t r,c;
	double v;

	while( ifile >> r >> c >> v ) {
		triplets.push_back( T(r-1,c-1,v) );
	}

	M->setFromTriplets( triplets.begin(), triplets.end() );
	ifile.close();
	return true;
}

bool readVector(const string& fname, std::vector<double>& v ) {
  
  ifstream ifile(fname);
  std::string line;
  for (size_t i = 0; i < 5; ++i) {
  	std::getline(ifile, line);
  }
  double ent;
  while ( ifile >> ent ) { v.push_back(ent); }
  ifile.close();
  return true;
}

    /**
    * Given the matrix A (m x n) and the right-hand-side vector b (m x 1) solve
    * for the solution vector x ( n x 1 ) satisfying
    * min_{x} || Ax - b ||
    * subject to x[i] >= 0 for all i
    */
    template< typename MatT, typename VecT >
    double* _solveNonNegativeLeastSquares( MatT& A, VecT& b ) {
        std::cerr << "in least squares problem\n";
        
        // Create a new TAUCS matrix structure
        taucs_ccs_matrix* mat;
        mat = new taucs_ccs_matrix;//(taucs_ccs_matrix*)malloc(sizeof(taucs_ccs_matrix));

        mat->n = A.cols();
        mat->m = A.rows();
        mat->flags = TAUCS_DOUBLE;
        
        std::cerr << "inited TAUCS matrix\n";

        // Compress this matrix so we can steal / share the data
        A.makeCompressed();

        std::cerr << "compressed matrix\n";
        
        mat->colptr = A.outerIndexPtr();
        mat->rowind = A.innerIndexPtr();
        mat->values.d = A.valuePtr();

        double residualNorm{0.0};

        std::cerr << "called NNLS solver\n";
        auto x = t_snnls_fallback( mat, &b[0], &residualNorm, 0.0, 0);
        char *errString;
        tsnnls_error(&errString);
        std::cerr << "NNLS solver returned\n";
        std::cerr << "ERROR STRING " << errString << "\n";

        if( x == NULL ) {
          std::cerr << "AHHH, x is STILL NULL\n";
          std::ofstream Afile("errorA.mtx");
          Afile << "SPARSE\n";
          Afile << A.rows() << '\t' << A.cols() << '\n' << A.nonZeros() << '\n';
          for (int k=0; k < A.outerSize(); ++k) {
            for (typename MatT::InnerIterator it(A,k); it; ++it) {
              Afile << it.row()+1 << '\t' << it.col()+1 << '\t' << it.value() << '\n';
            }
          }
          Afile.close();

          std::ofstream bfile("errorb.mtx");
          bfile << "# error right hand side\n";
          bfile << "# rows: " << b.size() << "\n";
          bfile << "# columns: 1\n";
          bfile << b.size() << '\t' << 1 << '\t' << b.size() << '\n';
          for (size_t i = 0; i < b.size(); ++i) {
              bfile << i << '\t' << 1 << '\t' << b[i] << '\n';
          }
          bfile.close();

        
          std::abort();
        }

        std::vector<double> myX( A.cols(), 0.0 );

        for( size_t i = 0; i < myX.size(); ++i ) { myX[i] = x[i]; std::cerr << "x[" << i << "]  = " << x[i] << "\n"; } 

        //free(x);
        delete mat;
      
        return x;//myX;
  
    }


int main(int argc, char* argv[]) {

	string mfile(argv[1]);
	string vfile(argv[2]);

	typedef Eigen::SparseMatrix<double> SpMat;
	SpMat* A;
	std::cerr << "reading matrix\n";
	readMatrix( mfile, A );
	std::cerr << "A is " << A->rows() << " x " << A->cols() << "\n";

	std::vector<double> b;
	std::cerr << "reading vector\n";
	readVector( vfile, b );

	auto x = _solveNonNegativeLeastSquares(*A,b);
	 
	if( x == NULL ) {
        	std::cerr << "DID NOT SOLVE\n";
        } else {
        	std::cerr << "SOLVED:\n";
        	for ( size_t i = 0; i < A->cols(); ++i ) {
        		std::cerr << " " << x[i];
        	}
        	std::cerr << "\n";
        }

	return 0;
}