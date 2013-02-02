#include <vector>
#include <fstream>
#include <memory>
#include <boost/range/irange.hpp>

#include "utils.hpp"
#include "matrix_tools.hpp"

int main(int argc, char* argv[]) {
	using std::vector;
	using std::string;
	using std::ifstream;

	ifstream ifile(argv[1]);
	auto t2g = utils::readTranscriptToGeneMap(ifile);
	std::cerr << "read " << t2g.numTranscripts() << " transcripts\n";
	std::cerr << "which mapped to " << t2g.numGenes() << " genes\n";
    ifile.close();


    /*
    *  A = [0 1 0 0 0 1 1]
	*      [0 0 1 1 0 1 1]
	*      [0 0 0 1 1 0 0]
	*  counts = [ 0 4 4 8 5 3 2 ]
	*/
    std::vector<std::vector<double>> A{{0,1,0,0,0,1,1},
                                       {0,0,1,1,0,1,1},
                                       {0,0,0,1,1,0,0}};
    std::vector<double> counts{0,4,4,8,5,3,2};

    std::unique_ptr<std::vector<std::vector<double>>> collapsedA{nullptr};
    std::unique_ptr<std::vector<double>> collapsedCounts{nullptr};
    matrix_tools::collapseIntoCategories(A, counts, collapsedA,collapsedCounts);

    for (size_t i = 0; i < collapsedA->size(); ++i ){
    	std::cerr << "[";
    	for ( size_t j = 0; j < (*collapsedA)[0].size(); ++j ) {
    		std::cerr << (*collapsedA)[i][j] << ", ";
    	}
    	std::cerr << "]\n";
    }
    std::cerr << "\n[";
    for ( auto c : *collapsedCounts ) { std::cerr << c << ", "; }
    std::cerr << "]\n";

    for ( auto i : boost::irange(0, 500) ) {
        std::cerr << i << ", ";
    }
}