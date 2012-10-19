#ifndef MATRIX_TOOLS_HPP
#define MATRIX_TOOLS_HPP

#include <unordered_set>

#include <Eigen/SparseCore>

namespace matrix_tools {

	template <typename MatT>
	std::vector<size_t> markDuplicateColumns( MatT& M ) {

		std::vector<size_t> removeCols;
		for ( size_t i = 0; i < M.cols(); ++i ) {
			Eigen::SparseVector<typename MatT::Scalar> colI(M.middleCols(i,1));
			auto normI = colI.dot(colI);
			auto stillValid = true;

			for ( size_t j = i+1; j < M.cols(); ++j ) {
				if ( stillValid ) {
					Eigen::SparseVector<typename MatT::Scalar> colJ(M.middleCols(j,1));
					auto normIJ = colI.dot(colJ);
					if ( normIJ > normI ) { 
						removeCols.push_back(i); 
						stillValid = false;
					}
				} 
			} // end j
		} // end i

		return removeCols;
	}

	/**
	* Project vector v onto vector u
	*/
	template <typename VecT>
	VecT projectOnto( VecT& v, VecT& u ) {
		return (v.dot(u) / u.dot(u) ) * u;
	}

	template <typename ColT>
	bool classicalGramSchmidt( std::vector<ColT>& basis, typename ColT::Scalar tol=1e-5 ) {

		// Orthogonalize the basis
		for ( size_t j = 0; j < basis.size(); ++j ) {
			// vj = aj
			ColT vj(basis[j]);
			for ( size_t i = 0; i < j; ++i ) {
				// rij = qi * aj
				//auto rij = basis[i].dot(basis[j]);
				// vj = vj - rii * qi
				vj = vj - projectOnto( basis[j], basis[i] );
				//(basis[i] * rij);
			}
			// rjj = ||vj||_{2}
			auto rjj = vj.dot(vj);
			// If this vector isn't orthogonal, then return false now
			if ( rjj <  tol ) { return false; }
			// qj = vj / rjj
			basis[j] = vj * (1.0 /  rjj);
		}

		return true;
	}

	template <typename ColT>
	bool addColumnToBasis( std::vector<ColT>& basis, ColT& c, typename ColT::Scalar tol=1e-5 ) {

		ColT vj(c);//basis[j]);
	    for ( size_t i = 0; i < basis.size(); ++i ) {
			// rij = qi * aj
	    	//auto rij = basis[i].dot(c);
			// vj = vj - rii * qi
	    	vj = vj - projectOnto( c, basis[i] );
	    	//vj = vj - (basis[i] * rij);
	    }
		// rjj = ||vj||_{2}
	    auto rjj = vj.dot(vj);

	    if ( rjj >= tol ) {
			// qj = vj / rjj
	    	basis.push_back( vj / rjj);
	    	return true;
	    } else {
	    	return false;
	    }
	}

	template <typename MatT>
	std::vector<size_t> markDependentColumns( MatT& M ) {

		typedef typename MatT::Index IndexT;
		typedef typename MatT::Scalar ScalarT;
		
		// Reference from MathOverflow:
		// http://mathoverflow.net/questions/109868/finding-linearly-independent-columns-of-a-large-sparse-rectangular-matrix
		
		// Structure that holds the column index and it's height
		struct ColHeightT {
			IndexT col;
			IndexT height;
		};

		// Compute the "height" of each column, where the "height" of a column is the
		// index of the row where the first non-zero entry occurs
		std::vector< ColHeightT > heights;
		heights.reserve(M.cols());

		for ( size_t i = 0; i < M.cols(); ++i ) {
			typename MatT::InnerIterator it(M,i);
			heights.push_back( {i, it.row()} );
		}

		// Sort the columns by height
		std::sort( heights.begin(), heights.end(), 
			[]( const ColHeightT& a, const ColHeightT& b) -> bool { return a.height < b.height; } );

		// We'll always keep the first column
		std::unordered_set<size_t> independentCols{ heights.front().col };
		IndexT currHeight{ heights.front().height };

		// Pick at least one column from each height class
		for ( auto& ch : heights ) { 
			if ( ch.height > currHeight ) { independentCols.insert(ch.col); currHeight = ch.height; }
		}

		std::vector< Eigen::SparseVector< ScalarT > > basis;
		for ( auto c : independentCols ) { 
			Eigen::SparseVector<ScalarT> col(M.middleCols(c,1)); 
			assert(col.size() == M.rows());
			basis.push_back( col );
		}

		// orthogonalize the current basis
		auto isIndependent = classicalGramSchmidt(basis);
	    if (!isIndependent) {
	    	std::cerr << "Something went horribly wrong; the original set of vectors was not linearly independent!\n";
	    	std::abort();
	    }

	    std::vector<size_t> removeCols;
	    // For each column that's not already in our basis
	    for ( size_t i = 0; i < M.cols(); ++i ) {
	    	if ( independentCols.find(i) == independentCols.end() ) {
	    		// Try to add it to our basis
	    		Eigen::SparseVector<ScalarT> c( M.middleCols(i,1) );
	    		assert(c.size() == M.rows());
	    		auto independent = addColumnToBasis(basis, c);
	    		// If we were able to add it, it's independent
	    		if ( independent ) { 
	    			independentCols.insert(i); 
	    		} else { // Otherwise it's dependent
	    			removeCols.push_back(i);
	    		}
	    	} // don't mess with already independent columsn
	    } // done loop over columns

	    return removeCols;
	}
}


#endif // MATRIX_TOOLS_HPP