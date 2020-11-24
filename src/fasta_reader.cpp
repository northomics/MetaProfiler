#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <Rcpp.h>
#include <chrono>
#include <random> 
#include <ctime>
#include <regex>
#include <algorithm>

// [[Rcpp::plugins(cpp11)]]


std::ifstream::pos_type filesize(const char* filename)
{
  std::ifstream input(filename, std::ifstream::ate | std::ifstream::binary);
  return input.tellg();
}
// [[Rcpp::export]]
Rcpp::DataFrame get_accessions(const std::vector<std::string> & files, const std::vector<std::string> & proteins)
{
  std::set<std::string> set;
  for(auto protein : proteins)
  {
    set.insert(protein);
  }
  time_t interval = time(nullptr);
  std::vector<std::string> accessions;
  std::vector<std::string> descriptors;
  std::vector<std::string> sequences;
  int fsize = 0;
  int bytes = 0;
  for (auto f : files)
  {
    fsize += filesize(f.c_str());
  }
  for (auto f : files)
  {
    std::ifstream file(f);
    std::string entry;
    getline(file, entry, '>');
    while(getline(file, entry, '>'))
    {
      if(set.find(entry.substr(0, entry.find_first_of(' '))) != set.end())
      {
        while (entry.find('\n') == std::string::npos) {
          std::string tmp = entry;
          getline(file, entry, '>');
          entry = tmp + '>' + entry;
        }
        accessions.push_back(entry.substr(0, entry.find_first_of(' ')));
        descriptors.push_back(entry.substr(entry.find_first_of(' ') + 1, entry.find_first_of('\n') - entry.find_first_of(' ') - 1));
        std::string seq = entry.substr(entry.find_first_of('\n') + 1, entry.size() - 1);
        seq.erase(std::remove(seq.begin(), seq.end(), '\n'), seq.end());
        sequences.push_back(seq);
      }
      bytes += entry.length();
      if (interval != time(nullptr))
      {
        interval = time(nullptr);
        Rcpp::Rcout << "reading fasta file(s) ... " << std::round(bytes * 100. / fsize) << "%  \r";
      }
    }
  }
  Rcpp::Rcout << "reading fasta file(s) ... 100%" << std::endl;
  return Rcpp::DataFrame::create(Rcpp::_["accession"] = accessions,
                                 Rcpp::_["description"] = descriptors,
                                 Rcpp::_["sequence"] = sequences,
                                 Rcpp::_["stringsAsFactors"] = false);
}


// [[Rcpp::export]]
Rcpp::DataFrame read_fasta(const std::vector<std::string> & files)
{
  time_t interval = time(nullptr);
  std::vector<std::string> accessions;
  std::vector<std::string> descriptors;
  std::vector<std::string> sequences;
  int fsize = 0;
  int bytes = 0;
  for (auto f : files)
  {
    fsize += filesize(f.c_str());
  }
  for (auto f : files)
  {
    std::ifstream file(f);
    std::string entry;
    getline(file, entry, '>');
    while(getline(file, entry, '>'))
    {
      while (entry.find('\n') == std::string::npos) {
        std::string tmp = entry;
        getline(file, entry, '>');
        entry = tmp + '>' + entry;
      }
      accessions.push_back(entry.substr(0, entry.find_first_of(' ')));
      descriptors.push_back(entry.substr(entry.find_first_of(' ') + 1, entry.find_first_of('\n') - entry.find_first_of(' ') - 1));
      std::string seq = entry.substr(entry.find_first_of('\n') + 1, entry.size() - 1);
      seq.erase(std::remove(seq.begin(), seq.end(), '\n'), seq.end());
      sequences.push_back(seq);
      bytes += entry.length();
      if (interval != time(nullptr))
      {
        interval = time(nullptr);
        Rcpp::Rcout << "reading fasta file(s) ... " << std::round(bytes * 100. / fsize) << "%  \r";
      }
    }
  }
  Rcpp::Rcout << "reading fasta file(s) ... 100%" << std::endl;
  return Rcpp::DataFrame::create(Rcpp::_["accession"] = accessions,
                                 Rcpp::_["description"] = descriptors,
                                 Rcpp::_["sequence"] = sequences,
                                 Rcpp::_["stringsAsFactors"] = false);
}


// [[Rcpp::export]]
void specific_fasta(const std::vector<std::string> & db, const std::vector<std::string> & proteins, std::string filename)
{
  std::set<std::string> set;
  for(auto protein : proteins)
  {
    set.insert(protein);
  }
  time_t interval = time(nullptr);
  int fsize = 0;
  int bytes = 0;
  for (auto f : db)
  {
    fsize += filesize(f.c_str());
  }
  std::ofstream ofile(filename.c_str());
  for (auto f : db)
  {
    std::ifstream file(f);
    std::string entry;
    getline(file, entry, '>');
    while(getline(file, entry, '>'))
    {
      if(set.find(entry.substr(0, entry.find_first_of(' '))) != set.end())
      {
        while (entry.find('\n') == std::string::npos) {
          std::string tmp = entry;
          getline(file, entry, '>');
          entry = tmp + '>' + entry;
        }
        ofile << '>' << entry;
      }
      bytes += entry.length();
      if (interval != time(nullptr))
      {
        interval = time(nullptr);
        Rcpp::Rcout << "reading fasta file(s) ... " << std::round(bytes * 100. / fsize) << "%  \r";
      }
    }
  }
  ofile.close();
  Rcpp::Rcout << "reading fasta file(s) ... 100%" << std::endl;
}

std::vector<std::string> split(const std::string & s, unsigned int missed_cleavage,
                               unsigned int min_length, unsigned int max_length)
{
  std::regex rgx("[RK]?[^P]\\w+?[KR](?!P)");
  std::sregex_token_iterator iter(s.begin(), s.end(), rgx);
  std::vector<std::string> res;
  std::sregex_token_iterator end;
  
  for ( ; iter != end; ++iter) {
    if((*iter).length() >= min_length && (*iter).length() <= max_length) {
      res.push_back(*iter);
    }
    auto old = iter;
    std::stringstream ss;
    ss << *iter;
    for(unsigned int i = 0; i < missed_cleavage; i++) {
      ++old;
      if(old != end) {
        ss << *old;
        if(ss.str().length() >= min_length && ss.str().length() <= max_length) {
          res.push_back(ss.str());
        }
      } else {
        break;
      }
    }
  }
  return res;
}



// [[Rcpp::export]]
std::set<std::string> trypsin_digestion(const std::vector<std::string> & files, unsigned int missed_cleavage,
                                        unsigned int min_length, unsigned int max_length) {
  time_t interval = time(nullptr);
  int fsize = 0;
  int bytes = 0;
  for (auto f : files)
  {
    fsize += filesize(f.c_str());
  }
  std::set<std::string> sequence;
  for (auto f : files)
  {
    std::ifstream file(f);
    std::string entry;
    getline(file, entry, '>');
    while(getline(file, entry, '>'))
    {
      while (entry.find('\n') == std::string::npos) {
        std::string tmp = entry;
        getline(file, entry, '>');
        entry = tmp + '>' + entry;
      }
      std::string seq = entry.substr(entry.find_first_of('\n') + 1, entry.size() - 1);
      seq.erase(std::remove(seq.begin(), seq.end(), '\n'), seq.end());
      std::vector<std::string> sseq = split(seq, missed_cleavage, min_length, max_length);
      for(auto s : sseq) {
        sequence.insert(s);
      }
      bytes += entry.length();
      if (interval != time(nullptr))
      {
        interval = time(nullptr);
        Rcpp::Rcout << "digesting sequences ... " << std::round(bytes * 100. / fsize) << "%  \r";
      }
    }
  }
  Rcpp::Rcout << "digesting sequences ... 100%" << std::endl;
  return sequence;
}

// [[Rcpp::export]]
void create_random_fasta(
    const std::vector<std::string> & files, std::string filename, const unsigned int db_size,
    const std::string add_with = "", const bool add_decoy = true,
    const std::string decoy_prefix = "DECOY_"
) {
  time_t interval = time(nullptr);
  
  
  
  
  unsigned int add_with_size = 0;
  bool contains_decoys = false;
  std::ofstream ofile(filename.c_str());
  std::set<std::string> sequences;
  if(add_with.length()) {
    std::ifstream file(add_with);
    std::string entry;
    getline(file, entry, '>');
    while(getline(file, entry, '>'))
    {
      while (entry.find('\n') == std::string::npos) {
        std::string tmp = entry;
        getline(file, entry, '>');
        entry = tmp + '>' + entry;
      }
      std::string header = entry.substr(0, entry.find_first_of('\n') + 1);
      std::string seq = entry.substr(entry.find_first_of('\n') + 1);
      seq.erase(std::remove(seq.begin(), seq.end(), '\n'), seq.end());
      if(entry.find(decoy_prefix) == std::string::npos) {
        add_with_size++;
        sequences.insert(seq);
        ofile << '>' << header << seq << '\n';
        if(add_decoy) {
          auto rev_seq = seq;
          std::reverse(rev_seq.begin(), rev_seq.end());
          ofile << ">DECOY_" << header << rev_seq << '\n';
        }
        
      } else {
        contains_decoys = true;

      }

    }
  }
  if(contains_decoys && add_decoy) {
    Rcpp::Rcout << "`add_with` already contains decoys. These will be added to the fasta file.\n";
  } else if(contains_decoys && !add_decoy) {
    Rcpp::Rcout << "`add_with` contains decoys. These will not be added to the fasta file.\n";
  }
  
  
  contains_decoys = false;
  int e = 0;
  std::vector<unsigned int> indices;
  for (auto f : files)
  {
    std::ifstream file(f);
    std::string entry;
    getline(file, entry, '>');
    while(getline(file, entry, '>'))
    {
      while (entry.find('\n') == std::string::npos) {
        std::string tmp = entry;
        getline(file, entry, '>');
        entry = tmp + '>' + entry;
      }
      if(entry.find(decoy_prefix) == std::string::npos) {
        indices.push_back(e++);
      } else {
        e++;
        contains_decoys = true;
      }
    }
  }
  
  if(contains_decoys && add_decoy) {
    Rcpp::Rcout << "The file(s) provided already contain decoys. These will be removed and new decoys will be generated.\n";
  }
  
  Rcpp::Rcout << indices.size() << " non-decoy entries in db(s)\n";
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::shuffle(indices.begin(), indices.end(), std::default_random_engine(seed));
  if(db_size < add_with_size) {
    Rcpp::stop("Number of entries in add_with is larger than the requested size.");
  }
  if(indices.size() < (db_size - add_with_size)) {
    Rcpp::stop("Requested size is too large for the number of entries in the given files.");
  }
  std::set<unsigned int> ran_indices(indices.begin(), indices.begin() + db_size - add_with_size);
  
  
  
  unsigned int ind = 0;
  for (auto f : files)
  {
    std::ifstream file(f);
    std::string entry;
    getline(file, entry, '>');
    while(getline(file, entry, '>'))
    {
      if(ran_indices.find(ind) == ran_indices.end()) {
        ind++;
        continue;
      }
      ind++;
      while (entry.find('\n') == std::string::npos) {
        std::string tmp = entry;
        getline(file, entry, '>');
        entry = tmp + '>' + entry;
      }
      std::string header = entry.substr(0, entry.find_first_of('\n') + 1);
      std::string seq = entry.substr(entry.find_first_of('\n') + 1);
      seq.erase(std::remove(seq.begin(), seq.end(), '\n'), seq.end());
      if(sequences.find(seq) == sequences.end()) {
        sequences.insert(seq);
        ofile << '>' << header << seq << '\n';
        if(add_decoy) {
          auto rev_seq = seq;
          std::reverse(rev_seq.begin(), rev_seq.end());
          ofile << ">DECOY_" << header << rev_seq << '\n';
        }
      } else {
        if((ind + 1) < (db_size - add_with_size)) {
          ran_indices.insert(ind + 1);
        }
      }
      if (interval != time(nullptr))
      {
        interval = time(nullptr);
        Rcpp::Rcout << "Creating fasta file ... " << sequences.size() * 100/db_size << "%  \r";
      }
      if(sequences.size() >= db_size) {
        break;
      }
    }
  }
  ofile.close();
  Rcpp::Rcout << "Creating fasta file ... 100%" << std::endl;
}

// [[Rcpp::export]]
int fasta_size(const std::vector<std::string> & files) {
  int num_seq = 0;
  for (auto f : files)
  {
    std::ifstream file(f);
    std::string entry;
    getline(file, entry, '>');
    while(getline(file, entry, '>'))
    {
      while (entry.find('\n') == std::string::npos) {
        std::string tmp = entry;
        getline(file, entry, '>');
        entry = tmp + '>' + entry;
      }
      ++num_seq;
    }
  }
  return num_seq;
}