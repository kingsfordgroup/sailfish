/* -*- Mode: C++; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/* vim: set ts=2 et sw=2 tw=80: */
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/filesystem.hpp>

#include <seqan/sequence.h>

#include "affymetrix_utils.hpp"

int main(int argc, char** argv) {

  using affymetrix::utils::AffyEntry;
  using std::vector;
  using std::string;

  std::string fname(argv[1]);

  auto pgmap = affymetrix::utils::ComputeReverseMap( fname );

  std::cerr << "size = " << pgmap.size() << "\n";


    std::string probeListFile(argv[2]);
    std::fstream ofile(argv[3], std::ios_base::out | std::ios_base::binary);
    auto probeFiles = affymetrix::utils::GetProbeFiles( probeListFile );

    vector<AffyEntry> ents;

    uint32_t probeLength = seqan::length(pgmap.begin()->second.probeSeq);
    uint64_t numRecords = 100;
    // Write a blank unit32_t and uint64_t to the file
    ofile.write( reinterpret_cast<char*>(&numRecords), sizeof(numRecords) );
    ofile.write( reinterpret_cast<char*>(&probeLength), sizeof(probeLength) );

    for ( auto& pf : probeFiles ) {
      std::cerr << "reading file " << pf << "\n";
      std::ifstream file(pf.native(), std::ios_base::in | std::ios_base::binary);
      try {

        boost::iostreams::filtering_istream in;
        in.push(boost::iostreams::gzip_decompressor());
        in.push(boost::iostreams::file_source(pf.native()));

        std::string str;
        std::getline(in, str);
        std::getline(in, str);

        while ( in.good() ) {

          std::string probeName;
          float expression, background;
          in >> probeName >> expression >> background;

          auto probeVal = pgmap.find( probeName );
          if ( probeVal != pgmap.end() ) {
            std::string probeSeq( seqan::String<char, seqan::CStyle>( pgmap[ probeName ].probeSeq ) );
            assert(probeSeq.size() == probeLength);
            ofile.write( probeSeq.c_str(), probeSeq.size()) ;
            ofile.write( reinterpret_cast<char*>(&expression), sizeof(expression));
            ofile.write( reinterpret_cast<char*>(&background), sizeof(background));
            ents.push_back( AffyEntry(probeName, expression, background) );
          }
        }

      } catch(const boost::iostreams::gzip_error& e) {
        std::cout << e.what() << '\n';
      }

      file.close();
    } // end for

    ofile.seekp(0);
    numRecords = ents.size();
    ofile.write( reinterpret_cast<char*>(&numRecords), sizeof(numRecords) );
    ofile.close();

    std::cout << "read " << ents.size() << " entries\n";
}
