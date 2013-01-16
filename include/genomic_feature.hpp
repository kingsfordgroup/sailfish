#ifndef GENOMIC_FEATURE_HPP
#define GENOMIC_FEATURE_HPP

#include <fstream>
#include <boost/tokenizer.hpp>

struct TranscriptGeneID {
  std::string transcript_id;
  std::string gene_id;
  // std::unordered_map<std::string, std::string> dynamic;
  bool parseAttribute( std::string& key, std::string& val ) {
    if ( key == "transcript_id" ) { transcript_id = val; return true; }
    if ( key == "gene_id" ) { gene_id = val; return true; }
    // dynamic[key] = val;
    return false;
  }

};

std::ostream& operator<< ( std::ostream& out, const TranscriptGeneID& ids ) {
  out << "transcript_id => " << ids.transcript_id << ", gene_id => " << ids.gene_id;
  return out;
}

template< typename StaticAttributes >
class GenomicFeature {
public:
   std::string seqID;
   std::string source;
   std::string type;
   int start, end;
   float score;
   char strand;
   char phase;
   StaticAttributes sattr;
   template <typename T>
   friend std::ostream& operator<< ( std::ostream& out, const GenomicFeature<T>& gf );
   template <typename T>
   friend std::istream& operator>> ( std::istream& in, const GenomicFeature<T>& gf );
};

template< typename StaticAttributes>
std::ostream& operator<< ( std::ostream& out, const GenomicFeature<StaticAttributes>& gf ) {
   out << "Feature: [" <<
          "id :" << gf.seqID << ", "  <<
          "source : " << gf.source << ", " <<
          "type : " << gf.type << ", " <<
          "start : " << gf.start << ", " <<
          "end : " << gf.end << ", " <<
          "score : " << gf.score << ", " <<
          "strand : " << gf.strand << ", " <<
          "phase : " << gf.phase << " ** ";
   out << gf.sattr << ", ";
   for ( auto& kv : gf.attributes ) {
      out << kv.first << " => " << kv.second << ", ";
   }
   return out;
}

template< typename StaticAttributes >
std::istream& operator>> ( std::istream& in, GenomicFeature<StaticAttributes>& gf ) {
  std::string line;
  std::getline(in, gf.seqID, '\t');
  std::getline(in, gf.source, '\t');
  std::getline(in, gf.type, '\t');
  
  std::getline(in, line, '\t');
  gf.start = atoi(line.c_str());

  std::getline(in, line, '\t');
  gf.end = atoi(line.c_str());

  std::getline(in, line, '\t');
  gf.score = atoi(line.c_str());

  std::getline(in, line, '\t');
  gf.strand = line[0];

  std::getline(in, line, '\t');
  gf.phase = line[0];

  
  std::getline(in, line);
  typedef boost::tokenizer<boost::char_separator<char> > 
    tokenizer;
  boost::char_separator<char> sep(";");  
  tokenizer tokens(line, sep);

  for ( auto tokIt : tokens ) {       
    // Currently, we'll handle the following separators
    // '\s+'
    // '\s*=\s*'
    tokIt = tokIt.substr(tokIt.find_first_not_of(' '));
    auto kvsepStart = tokIt.find('=');
    
    // If we reached the end of the key, value token, then the string must have been
    // separated by some set of spaces, and NO '='.  If this is the case, find the 'spaces' so that
    // we can split on it.
    if ( kvsepStart == tokIt.npos ) {
      kvsepStart = tokIt.find(' ');
    } 
    
    auto key = tokIt.substr(0, kvsepStart);
    key = key.substr(0, key.find(' '));

    auto kvsepStop = 1 + kvsepStart + tokIt.substr(kvsepStart+1).find_first_not_of(' ');
    auto val = (tokIt[kvsepStop] == '"') ? tokIt.substr(kvsepStop+1, (tokIt.length() - (kvsepStop+2))) : tokIt.substr(kvsepStop, (tokIt.length() - (kvsepStop+1)));
    gf.sattr.parseAttribute( key, val );
  }

  return in;
}

#endif // GENOMIC_FEATURE_HPP