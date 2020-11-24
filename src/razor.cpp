#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <vector>
#include <string>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List razor(std::vector<std::vector<std::string>> x, std::vector<std::string> id, std::vector<unsigned int> unique, const bool verbose = true)
{
  if(x.size() != id.size() || x.size() != unique.size() || unique.size() != id.size()) {
    Rcpp::stop("x, id, and unique are not the same size.");
  }
  
  std::vector<unsigned int> sorted_indices(x.size());
  for(unsigned int i=0;i<x.size();i++)
    sorted_indices[i] = i;
  std::sort(sorted_indices.begin(), sorted_indices.end(),
            [&](const unsigned int& l, const unsigned int& r) {
              if(x[l].size() == x[r].size()) {
                return unique[l] > unique[r];
              }
              return x[l].size() > x[r].size();
            }
  );
  
  std::vector<std::vector<std::string>> x_sorted(x.size());
  std::vector<std::string> id_sorted(id.size());
  std::vector<unsigned int> unique_sorted(unique.size());
  for(unsigned int i=0;i<sorted_indices.size();i++)
  {
    x_sorted[i] = x[sorted_indices[i]];
    id_sorted[i] = id[sorted_indices[i]];
    unique_sorted[i] = unique[sorted_indices[i]];
  }
  
  std::vector<std::vector<std::string>> group(1);
  group[0].push_back(id_sorted[0]);
  std::vector<std::vector<std::string>> cmp;
  cmp.push_back(x_sorted[0]);
  std::vector<std::vector<unsigned int>> count(1);
  count[0].push_back(x_sorted[0].size());
  std::vector<std::vector<unsigned int>> unique_(1);
  unique_[0].push_back(unique_sorted[0]);
  time_t last_invoke_ = time(nullptr);
  for(unsigned int i = 1; i < x_sorted.size(); i++)
  {
    for(unsigned int j = 0; j < cmp.size(); j++)
    {
      for(std::vector<std::string>::iterator it = x_sorted[i].begin(); it != x_sorted[i].end(); it++)
      {
        auto found = find(cmp[j].begin(), cmp[j].end(), *it);
        if(found == cmp[j].end())
        {
          goto ctn1;
        }
      }
      group[j].push_back(id_sorted[i]);
      count[j].push_back(x_sorted[i].size());
      unique_[j].push_back(unique_sorted[i]);
      goto ctn2;
      ctn1:;
    }
    cmp.push_back(x_sorted[i]);
    group.push_back(std::vector<std::string>(1, id_sorted[i]));
    count.push_back(std::vector<unsigned int>(1, x_sorted[i].size()));
    unique_.push_back(std::vector<unsigned int>(1, unique_sorted[i]));
    ctn2:;
    if (verbose && last_invoke_ != time(nullptr))
    {
      last_invoke_ = time(nullptr);
      Rcpp::Rcout << round((double)(i + 2)/(double)x_sorted.size() * 100) << "%                     \r";
    }
  }
  if (verbose)
  {
    Rcpp::Rcout << "100%                     \n";
  }
  return Rcpp::List::create(Rcpp::_["id"] = group, Rcpp::_["count"] = count, Rcpp::_["unique"] = unique_, Rcpp::_["x"] = cmp);
}
