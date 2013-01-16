#ifndef MATRIX_TOOLS_HPP
#define MATRIX_TOOLS_HPP

#include <unordered_set>
#include <unordered_map>
#include <map>
#include <memory>

#include <Eigen/SparseCore>
#include <boost/dynamic_bitset.hpp>

namespace matrix_tools {

/**
*  This function collapses the rate matrix "A" and the corresponding vector of read
*  counts "counts" into a set of unique categories as defined in the paper 
*  (Salzman, Jiang and Wong 2011).

*  Example 
*  A = [0 1 0 0 0 1 1]
*      [0 0 1 1 0 1 1]
*      [0 0 0 1 1 0 0]
*  counts = [ 0 4 4 8 5 3 2 ]
* 
*  might (b/c output column order is undefined) be collapsed to:
*  A = [1 0 0 0 2]
*      [0 1 1 0 2]
*      [0 0 1 1 0]
*  counts = [4 4 8 5 5]
*/
void collapseIntoCategories(
	std::vector<std::vector<double>> &A,
	std::vector<double> &counts,
	std::unique_ptr<std::vector<std::vector<double>>>& collapsedA,
	std::unique_ptr<std::vector<double>>& collapsedCounts
) {

    using boost::dynamic_bitset;
    typedef size_t Count;
    typedef size_t Index;

    size_t numRows = A.size();
    size_t numCols = A[0].size();
    std::map<dynamic_bitset<>, std::vector<Index> > categoryCounts;

    for ( size_t i = 0; i < numCols; ++i ) {
        dynamic_bitset<> category(numRows);

        for ( size_t j = 0; j < numRows; ++j ) {
            category[j] = ( A[j][i] > 0.0 ) ? 1 : 0;
        }
        // The trivial category contains all 0s; we simply omit it
        bool trivial = !category.any();
        if ( !trivial ) {
            categoryCounts[ category ].push_back(i);
        }
    }

    size_t numCategories = categoryCounts.size();
    std::unique_ptr<std::vector<std::vector<double>>> captr(new std::vector<std::vector<double>>(
    	numRows, std::vector<double>( numCategories, 0.0 ) ));
    collapsedA = std::move(captr);

    std::unique_ptr<std::vector<double>> ccptr( new std::vector<double>(numCategories, 0.0) );
    collapsedCounts = std::move(ccptr);

    size_t j=0;
    for ( auto catIt = categoryCounts.begin(); catIt != categoryCounts.end(); ++catIt, ++j ) {
        auto category = catIt->first;
        auto origCols = catIt->second;
        for ( auto origColID : origCols ) {
            (*collapsedCounts)[j] += counts[origColID];
        }
        for ( size_t i = 0; i < numRows; ++i ) {
            (*collapsedA)[i][j] = origCols.size() * category[i];
        }
    }

}



template <typename MatT>
std::vector<size_t> markDuplicateColumns( MatT &M ) {

    std::vector<size_t> removeCols;
    for ( size_t i = 0; i < M.cols(); ++i ) {
        Eigen::SparseVector<typename MatT::Scalar> colI(M.middleCols(i, 1));
        auto normI = colI.dot(colI);
        auto stillValid = true;

        for ( size_t j = i + 1; j < M.cols(); ++j ) {
            if ( stillValid ) {
                Eigen::SparseVector<typename MatT::Scalar> colJ(M.middleCols(j, 1));
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
VecT projectOnto( VecT &v, VecT &u ) {
    return (v.dot(u) / u.dot(u) ) * u;
}

template <typename ColT>
bool classicalGramSchmidt( std::vector<ColT> &basis, typename ColT::Scalar tol = 1e-10 ) {

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
        if ( rjj <  tol ) {
            return false;
        }
        // qj = vj / rjj
        basis[j] = vj * (1.0 /  rjj);
    }

    return true;
}

template <typename ColT>
bool addColumnToBasis( std::vector<ColT> &basis, ColT &c, typename ColT::Scalar tol = 1e-10 ) {

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
std::vector<size_t> markDependentColumns( MatT &M ) {

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
        typename MatT::InnerIterator it(M, i);
        heights.push_back( {i, it.row()} );
    }

    // Sort the columns by height
    std::sort( heights.begin(), heights.end(),
               []( const ColHeightT & a, const ColHeightT & b) -> bool { return a.height < b.height; } );

    // We'll always keep the first column
    std::unordered_set<size_t> independentCols { heights.front().col };
    IndexT currHeight { heights.front().height };

    // Pick at least one column from each height class
    for ( auto & ch : heights ) {
        if ( ch.height > currHeight ) {
            independentCols.insert(ch.col);
            currHeight = ch.height;
        }
    }

    std::vector< Eigen::SparseVector< ScalarT > > basis;
    for ( auto c : independentCols ) {
        Eigen::SparseVector<ScalarT> col(M.middleCols(c, 1));
        assert(col.size() == M.rows());
        basis.push_back( col );
    }

    // orthogonalize the current basis
    auto isIndependent = classicalGramSchmidt(basis);
    if (!isIndependent) {
        std::cerr << "Something went horribly wrong; the original set of vectors was not linearly independent!\n";
        std::cerr << "Doing the safe thing (start with just one col and build up\n";
        independentCols = {0};
        Eigen::SparseVector<ScalarT> col(M.middleCols(0, 1));
        basis = { col };
        classicalGramSchmidt(basis);
        //std::abort();
    }

    std::vector<size_t> removeCols;
    // For each column that's not already in our basis
    for ( size_t i = 0; i < M.cols(); ++i ) {
        if ( independentCols.find(i) == independentCols.end() ) {
            // Try to add it to our basis
            Eigen::SparseVector<ScalarT> c( M.middleCols(i, 1) );
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