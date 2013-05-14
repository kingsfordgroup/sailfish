/* -*- Mode: C++; tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/* vim: set ts=2 et sw=2 tw=80: */

#include "affymetrix_utils.hpp"
#include <seqan/seq_io.h>

namespace affymetrix {

  namespace utils {

    using std::string;
    using std::vector;

    std::ostream& operator<< (std::ostream& os, const ProbeTarget& p) {
      os << "probe target: name [" << p.geneName << "], seq [" << p.probeSeq << "]";
      return os;
    }

    ProbeTarget::ProbeTarget( const string& _geneName, const seqan::DnaString& _probeSeq ) : geneName(_geneName), probeSeq(_probeSeq){}
    ProbeTarget::ProbeTarget() : geneName(""), probeSeq(""){}
   /**
    * We want to be able to hash
    */

    /**
     * Given a file containing the probe information, build a reverse mapping from
     * probe name to ENSEMBL gene name, and probe sequence.  Any probe mapping to > 1 gene will be discarded.
     */
    StringMapT ComputeReverseMap( const string& fname ) {
      StringMapT pgmap;
      StringSetT duplicates;

      std::ifstream in( fname.c_str() );

      auto contains = [&]( const string& k ) -> bool { return !(pgmap.find(k) == pgmap.end());};

      // Throw away the first line
      string line;
      std::getline(in, line );

      string probe, gene;

      boost::char_separator<char> sep(",", "", boost::keep_empty_tokens);

      for( ; std::getline(in, line); ) {
        boost::tokenizer<boost::char_separator<char>> tok(line, sep);
        auto tokit = tok.begin();
        gene =  *tokit; ++tokit; probe = *tokit;
        for( size_t i = 0; i < 8; ++i ) { ++tokit; }
        seqan::DnaString seq = *tokit;

        if ( !contains(gene) ) {
          pgmap[probe] = ProbeTarget{gene, seq};
        } else {
          duplicates.insert( probe );
        }
      }

      for ( auto& dup : duplicates ) { pgmap.erase(dup); }

      return pgmap;
    }


    /**
     * Given a file containing a list of strings, return a vector of corresponding
     * paths.  This is the list of probe files to be read.
     */
    vector< boost::filesystem::path > GetProbeFiles( const string& fname ) {

      vector< boost::filesystem::path > probeFiles;

      boost::filesystem::path dirPath(fname);
      dirPath = dirPath.parent_path();

      std::ifstream in(fname.c_str());
      for ( std::string str; std::getline(in,str); ) {
        boost::filesystem::path p(dirPath);
        p /= boost::filesystem::path(str);
        assert(boost::filesystem::exists(p));
        probeFiles.push_back(p);
      }

      return probeFiles;
    }

    AffyEntry::AffyEntry( const string& _probeName, float _expression, float _background ) :
      probeName(_probeName), expression(_expression), background(_background) {}
  }

}
