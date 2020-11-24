#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <map>
#include <set>
#include <ctime>
#include <string>
#include <vector>
#include <iostream>
#include <mutex>
#include <cmath>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]


inline bool cols_equal(const arma::uvec & lhs, const arma::uvec & rhs) {
  return arma::all(lhs == rhs);
}

arma::uvec unique_cols(const arma::umat& x, arma::uvec & indices) {
  unsigned int count = 1, i = 1, j = 1, nr = x.n_rows, nc = x.n_cols;
  arma::umat result(nr, nc);
  result.col(0) = x.col(0);
  
  indices(0) = 0;
  
  for ( ; i < nc; i++) {
    bool matched = false;
    if (cols_equal(x.col(i), result.col(0))) continue;
    
    for (j = i + 1; j < nc; j++) {
      if (cols_equal(x.col(i), x.col(j))) {
        matched = true;
        break;
      }
    }
    
    if (!matched) {
      indices(count) = i;
      result.col(count++) = x.col(i);
    }
  }
  return indices.rows(0, count - 1);
}

void unique_cols(arma::umat& x, arma::uvec & indices) {
  unsigned int count = 1, i = 1, j = 1, nr = x.n_rows, nc = x.n_cols;
  arma::umat result(nr, nc);
  result.col(0) = x.col(0);
  
  indices(0) = 0;
  
  for ( ; i < nc; i++) {
    bool matched = false;
    if (cols_equal(x.col(i), result.col(0))) continue;
    
    for (j = i + 1; j < nc; j++) {
      if (cols_equal(x.col(i), x.col(j))) {
        matched = true;
        break;
      }
    }
    
    if (!matched) {
      indices(count) = i;
      result.col(count++) = x.col(i);
    }
  }
  x = x.cols(0, count - 1);
  for(unsigned int i = 0; i < count; i++) {
    x.col(i) = result.col(i);
  }
}

void unique_cols(arma::umat& x) {
  unsigned int count = 1, i = 1, j = 1, nr = x.n_rows, nc = x.n_cols;
  arma::umat result(nr, nc);
  result.col(0) = x.col(0);
  for ( ; i < nc; i++) {
    bool matched = false;
    if (cols_equal(x.col(i), result.col(0))) continue;
    
    for (j = i + 1; j < nc; j++) {
      if (cols_equal(x.col(i), x.col(j))) {
        matched = true;
        break;
      }
    }
    
    if (!matched) {
      result.col(count++) = x.col(i);
    }
  }
  x = x.cols(0, count - 1);
  for(unsigned int i = 0; i < count; i++) {
    x.col(i) = result.col(i);
  }
}

int binary_search(arma::umat & x, unsigned int size, int low, int high) 
{ 
  if (high <= low)
    return (size <= arma::accu(x.col(low))) ? (low + 1) : low; 
  
  int mid = (low + high)/2;
  unsigned int mid_val = arma::accu(x.col(mid));
  if(size == mid_val)
  {
    return mid+1; 
  }
  
  if(size <= mid_val) 
    return binary_search(x, size, mid+1, high); 
  return binary_search(x, size, low, mid-1); 
} 

// Function to sort an array a[] of size 'n' 
void sort_mat(arma::umat & x, arma::uvec & indices, unsigned int start) 
{ 
  for (unsigned int i = start + 1; i < x.n_cols; ++i) 
  { 
    int j = i - 1; 
    unsigned int current_index = indices(i); 
    arma::uvec current_col = x.col(i); 
    // find location where selected sould be inseretd 
    int index = binary_search(x, arma::accu(current_col), start, j);
    // Move all elements after location to create space 
    while (j >= index) 
    {
      indices(j+1) = indices(j);
      x.col(j+1) = x.col(j); 
      j--; 
    } 
    indices(index) = current_index; 
    x.col(index) = current_col; 
  }
} 

void sort_mat(arma::umat & x, unsigned int start) 
{ 
  for (unsigned int i = start + 1; i < x.n_cols; ++i) 
  { 
    int j = i - 1; 
    arma::uvec current_col = x.col(i); 
    // find location where selected sould be inseretd 
    int index = binary_search(x, arma::accu(current_col), start, j);
    // Move all elements after location to create space 
    while (j >= index) 
    {
      x.col(j+1) = x.col(j); 
      j--; 
    } 
    x.col(index) = current_col; 
  }
} 


bool has_overlap(arma::umat & x) {
  arma::uvec check = arma::sum(x, 1);
  for(auto c : check) {
    if(c > 1) return true;
  }
  return false;
}
// 
// std::vector<std::vector<int>> permutations(std::vector<int> & x) 
// { 
//   std::vector<std::vector<int>> indices;
//   indices.push_back(x);
//   while (std::next_permutation(x.begin(), x.end())) {
//      indices.push_back(x);
//   }
//   return indices;
// } 
// 
// void get_indices(std::vector<std::vector<int>> & final,
//                  std::vector<int> & current,
//                  std::vector<std::vector<std::vector<int>>>::const_iterator first,
//                  std::vector<std::vector<std::vector<int>>>::const_iterator last) {
//   if (first == last) {
//     final.push_back(current);
//     return;
//   }
//   const std::vector<std::vector<int>> & tmp = *first;
//   for (std::vector<std::vector<int>>::const_iterator it = tmp.begin(); it != tmp.end(); it++) {
//     current.insert(current.end(), (*it).begin(), (*it).end());
//     get_indices(final, current, first+1, last);
//     current.resize(current.size() - (*it).size());
//   }
// }

arma::umat remove_empty_columns(arma::umat x) {
  int count = 0;
  arma::urowvec col_sums = arma::sum(x, 0);
  for(unsigned int i = 0; i < col_sums.n_cols; i++) {
    if(col_sums(i) == 0) {
      break;
    }
    count = i;
  }
  return x.cols(0, count);
}
// 
// void combn_(std::vector<std::vector<unsigned int>> & final, 
//             std::vector<unsigned int> & current, unsigned int n, unsigned int left, unsigned int k) 
// { 
//   if (k == 0) { 
//     final.push_back(current); 
//     return; 
//   } 
//   
//   for (unsigned int i = left; i <= n; ++i) 
//   { 
//     current.push_back(i); 
//     combn_(final, current, n, i + 1, k - 1); 
//     current.pop_back(); 
//   } 
// }

// std::vector<std::vector<unsigned int>> combn(unsigned int n, unsigned int k) 
// { 
//   std::vector<std::vector<unsigned int>> final; 
//   std::vector<unsigned int> current; 
//   combn_(final, current, n - 1, 0, k); 
//   return final; 
// } 

arma::mat relative_distance(const arma::colvec & x, const double y)
{
  return arma::sum(arma::abs(x - y)/(arma::abs(x + y)/2), 1)*100;
}

arma::mat absolute_distance(const arma::colvec & x, const double y)
{
  return arma::sum(arma::abs(x - y), 1);
}

arma::mat bray_distance(const arma::colvec & x, const double y)
{
  return arma::sum(arma::abs(x - y),1)/arma::sum(x + y, 1);
}

arma::mat jaccard_index(const arma::colvec & x, const double y)
{
  arma::mat B = bray_distance(x,y);
  return 2 * B/(1 + B);
}

double euclidean_distance(const arma::rowvec & x, const arma::rowvec & y)
{
  return std::sqrt(arma::accu(arma::pow(x - y, 2)));
}

double bray_distance(const arma::rowvec & x, const arma::rowvec & y)
{
  return arma::accu(arma::abs(x - y))/arma::accu(x + y);
}

double jaccard_index(const arma::rowvec & x, const arma::rowvec & y)
{
  double B = bray_distance(x,y);
  return 2 * B/(1 + B);
}

double correlation(const arma::rowvec & x, const arma::rowvec & y)
{
  arma::rowvec rhs = x - arma::mean(x);
  double rhs_squared_sum = arma::accu(arma::pow(rhs, 2));
  arma::rowvec lhs = y - arma::mean(y);
  double lhs_squared_sum = arma::accu(arma::pow(lhs, 2)); 
  return arma::accu(lhs % rhs)/std::sqrt(lhs_squared_sum * rhs_squared_sum);
}


double relative_distance(const arma::rowvec & x, const arma::rowvec & y)
{
  return arma::mean(arma::abs(x - y)/((x + y)/2));
}

void makeClusterMap( std::vector<arma::uvec> & x, std::vector<arma::uvec> & y)
{
  arma::ivec nums(x[0].size());
  for(unsigned int i = 1; i <= nums.size(); i++) {
    nums(i - 1) = i;
  }
  int max = nums.size() * (nums.size() + 1) / 2;
  y.push_back(x[0]);
  x.erase(x.begin());
  while(x.size() > 0)
  {
    for(unsigned int i = 0; i < x.size(); i++)
    {
      for(unsigned int j = 0; j < x[i].size(); j++)
      {
        if(x[i](j) && y[y.size() - 1](j))
        {
          x[i](j) = false;
        }
      }
    }
    std::sort(x.begin(), x.end(), [&](const arma::uvec & a, const arma::uvec & b) -> bool 
    { 
      int lhs = max * arma::accu(a) + arma::accu(a % nums);
      int rhs = max * arma::accu(b) + arma::accu(b % nums);
      return lhs > rhs;
    }) ;
    x.erase(std::unique(x.begin(), x.end(), [](const arma::uvec & a, const arma::uvec & b) -> bool
    {
      return arma::all(a == b);
    }), x.end());

    if(!arma::any(x[0]))
    {
      return;
    }
    y.push_back(x[0]);
    x.erase(x.begin());
  }
}

void dist_rowwise(const arma::mat & y,
                  std::vector<arma::uvec> & clusterMap,
                  const arma::mat & radius,
                  const Rcpp::StringVector & method,
                  const Rcpp::IntegerVector & start,
                  const Rcpp::IntegerVector & end)
{
  arma::uvec dok(y.n_rows);
  for(unsigned int i = 0;  i < y.n_rows; i++)
  {
    arma::mat d(radius.n_rows, radius.n_cols);
    for(unsigned int j = 0;  j < y.n_rows; j++)
    {
      if(j == i)
      {
        dok[j] = true;
        continue;
      }
      for (int c = 0; c < method.size(); c++)
      {
        if (method[c] == "correlation")
        {
          d(j,c) = 1 - correlation(y.row(j).cols(start[c], end[c]), y.row(i).cols(start[c], end[c]));
        } 
        else if (method[c] == "bray")
        {
          d(j,c) = bray_distance(y.row(j).cols(start[c], end[c]), y.row(i).cols(start[c], end[c]));
        }
        else if (method[c] == "jaccard")
        {
          d(j,c) = jaccard_index(y.row(j).cols(start[c], end[c]), y.row(i).cols(start[c], end[c]));
        }
        else if (method[c] == "relative")
        {
          d(j,c) = relative_distance(y.row(j).cols(start[c], end[c]), y.row(i).cols(start[c], end[c]));
        }
        else if(method[c] == "euclidean")
        {
          d(j,c) = euclidean_distance(y.row(j).cols(start[c], end[c]), y.row(i).cols(start[c], end[c]));
        }
      }
    }
    arma::mat cmp(y.n_rows, 1);
    cmp.fill(y.n_cols);
    auto z = d <= radius;
    auto w = arma::sum(z, 1);
    arma::uvec dok = w >= cmp;
    clusterMap.push_back(dok);
  }
}

int fact(int n)    
{    
  if(n<0)    
    return(-1); /*Wrong value*/      
  if(n==0)    
    return(1);  /*Terminating condition*/    
  else    
  {    
    return(n*fact(n-1));        
  }    
}

arma::umat greater_than(arma::mat a, arma::mat b, float epsilon)
{
  if((a.n_rows != b.n_rows) || (a.n_cols != b.n_cols)) {
    Rcpp::stop("a and b must have equal dimensions for greater than comparison.");
  }
  arma::umat res(a.n_rows, a.n_cols);
  for(unsigned int r = 0; r < a.n_rows; r++) {
    for(unsigned int c = 0; c < a.n_cols; c++) {
      double lhs = a(r,c);
      double rhs = b(r,c);
      res(r, c) = (lhs - rhs) > (std::max(std::fabs(lhs), std::fabs(rhs)) * epsilon);
    }
  }
  return res;
}

// void keep_valid_combn(std::vector<std::vector<unsigned int>> & combn_vec, const arma::umat & clusts) {
//   std::vector<std::vector<unsigned int>> tmp;
//   for(unsigned int i = 0; i < combn_vec.size(); i++) {
//     arma::umat current(clusts.n_rows, combn_vec[i].size());
//     for(unsigned int j = 0; j < combn_vec[i].size(); j++) {
//       current.col(j) = clusts.col(combn_vec[i][j]);
//     }
//     arma::uvec row_sums = arma::sum(current, 1);
//     bool valid = true;
//     for(unsigned int i = 0; i < row_sums.n_rows; i++) {
//       if(row_sums(i) == 0) {
//         valid = false;
//         break;
//       }
//     }
//     if(valid) {
//       tmp.push_back(combn_vec[i]);
//     }
//   }
//   combn_vec = tmp;
// }

// arma::umat get_clust() {
//   arma::umat tmp;
//   for(unsigned int i = 0; i < best_indices.size(); ++i) {
//     tmp.col(i) = clusts.col(best_indices[i]);
//   }
//   return tmp;
// }

// [[Rcpp::export]]
Rcpp::List qtclust_c (arma::mat & m,
                      const int n_groups,
                      arma::uvec & id,
                      arma::mat & groups,
                      arma::vec & radius, 
                      Rcpp::StringVector& method,
                      Rcpp::IntegerVector & start,
                      Rcpp::IntegerVector & end,
                      const bool merge_overlaps,
                      const bool element_wise,
                      const bool verbose)
{
  time_t last_invoke_ = time(nullptr);
  std::vector<arma::ivec> clusters(n_groups); // container for cluster id
  std::vector<arma::mat> centers(n_groups); // container for the means in each cluster
  std::vector<arma::mat> errors(n_groups); // container for the sd of the means in each cluster
  // std::vector<arma::rowvec> quality(n_groups); // container for silhouette measure
  unsigned int rows_done = 0;
  // unsigned int m_n_cols = m.n_cols;
  unsigned int gr_n_cols = groups.n_cols;
  int nb_clust = 0;
  
  int k = 1;
  
  unsigned int index = 0;
  while(rows_done < m.n_rows) {
    int prev_id = id(rows_done);
    int current_id = prev_id;
    unsigned int current_row = rows_done;
    while((prev_id == current_id) && (current_row < (m.n_rows - 1))) {
      current_id = id(++current_row);
    }
    if(current_row == (m.n_rows - 1)) {
      current_row++;
    }
    arma::mat y = m.rows(rows_done, current_row - 1);
    arma::mat gr = groups.rows(rows_done, current_row - 1);
    rows_done = current_row;
    if(y.n_rows == 1)
    {
      clusters[index] = k;
      arma::mat tmp(gr.n_rows, gr_n_cols + 1);
      tmp.cols(0, gr_n_cols - 1) = gr;
      tmp.col(gr_n_cols) = 1;
      centers[index] = tmp;
      errors[index] = arma::mat(1, gr_n_cols, arma::fill::zeros);
      // quality[index] = arma::rowvec(gr_n_cols, arma::fill::zeros);
      
      k++;
      nb_clust++;
      index++;
      
      continue;
    }
    std::vector<arma::mat> d_vec(y.n_rows);
    /// std::vector<arma::vec> d_max_vec(y.n_rows);
    if (element_wise) // distance measurements is done element-wise
    {
      
      // test each row as a potential cluster center and keep the one that cluster the most rows
      for(unsigned int i = 0;  i < y.n_rows; i++)
      {
        arma::mat d(y.n_rows, y.n_cols);
        for(unsigned int j = 0;  j < y.n_cols; j++) // distance measured element-wise
        {
          if(method[j] == "absolute")
          {
            d.col(j) = absolute_distance(y.col(j), y(i,j));
          }
          else if(method[j] == "relative")
          {
            d.col(j) = relative_distance(y.col(j), y(i,j));
          }
          else if(method[j] == "bray")
          {
            d.col(j) = bray_distance(y.col(j), y(i,j));
          }
          else if(method[j] == "relative")
          {
            d.col(j) = jaccard_index(y.col(j), y(i,j));
          }
          else
          {
            Rcpp::stop("not a valid distance method");
          }
        }
        
        // arma::vec res(m_n_cols);
        // for (unsigned int k = 0; k < d.n_cols; ++k) {
        //   res(k) = d.col(k).max();
        // }

        d_vec[i] = d;
        // d_max_vec[i] = res;
      }
     
      arma::umat clusts(y.n_rows, y.n_rows);
      arma::mat radius_mat(y.n_rows, y.n_cols);
      for(unsigned int i = 0;  i < y.n_rows; i++) {
        radius_mat.row(i) = radius.t();
      }
      int i = 0;
      for(auto d : d_vec) {
        auto z = greater_than(d, radius_mat, 1e-6);
        auto w = arma::sum(z, 1);
        arma::uvec dok = w == 0;
        clusts.col(i++) = dok;
      }

      // Rcpp::Rcout << clusts;
     
      arma::uvec indices(clusts.n_cols);
      unique_cols(clusts);
      // Rcpp::Rcout << clusts;
      sort_mat(clusts, 0);
      // Rcpp::Rcout << clusts;
      arma::umat best_clust = clusts;
      unsigned int best_size = clusts.n_cols;
      if(has_overlap(clusts)) {
        // Rcpp::Rcout << clusts;
        // Rcpp::Rcout << std::endl;
        // 
        // arma::mat distance_matrix(best_clust.n_cols, best_clust.n_cols);
        // for(unsigned int j = 0; j < best_clust.n_cols; j++) {
        //   distance_matrix(j, j) = 0.;
        //   for(unsigned int k = j + 1; k < best_clust.n_cols; k++) {
        //     arma::rowvec num = d_vec[indices(j)].row(indices(k));
        //     arma::rowvec denom(m.n_cols);
        //     for(unsigned int c = 0; c < m_n_cols; c++) {
        //       denom(c) = std::max(d_max_vec[indices(j)](c), d_max_vec[indices(k)](c));
        //     }
        //     arma::rowvec row_val = denom/num;
        //     double val = arma::accu(row_val)/m_n_cols;
        //     distance_matrix(j, k) = val;
        //     distance_matrix(k, j) = val;
        //   }
        // }
        // 
        // double best_result = -1.;
        // std::vector<unsigned int> best_indices;
        // 
        // for(unsigned int n_groups = 1; n_groups < clusts.n_cols; n_groups++) {
        //   std::vector<std::vector<unsigned int>> combn_vec = combn(clusts.n_cols, n_groups);
        //   keep_valid_combn(combn_vec, clusts);
        //   if(combn_vec.size() == 0) {
        //     continue;
        //   }
        //   for(auto i : combn_vec) {
        //     double sum = 0;
        //     for(unsigned int j = 0; j < i.size(); j++) {
        //       for(unsigned int k = j + 1; k < i.size(); k++) {
        //         sum += distance_matrix(j, k);
        //       }
        //     }
        //     if((best_result < 0) || (sum < best_result)) {
        //       best_indices = i;
        //       best_result = sum;
        //     }
        //   }
        // }
        // Rcpp::Rcout <<  "best result is " << best_result << std::endl;
        // 
        // for(unsigned int i = 0; i < best_indices.size(); ++i) {
        //   best_clust.col(i) = clusts.col(best_indices[i]);
        // }
        // best_clust = best_clust.cols(0, best_indices.size() - 1);
        // sort_mat(best_clust, 0);
        
        if(merge_overlaps)
        {
          for(unsigned int col1 = 0; col1 < best_clust.n_cols; col1++)
          {
            for(unsigned int row1 = 0; row1 < best_clust.n_rows; row1++)
            {
              if (best_clust(row1,col1))
              {
                for(unsigned int col2 = col1 + 1; col2 < best_clust.n_cols; col2++)
                {
                  if(best_clust(row1,col2))
                  {
                    for(unsigned int row2 = 0; row2 < best_clust.n_rows; row2++)
                    {
                      if(best_clust(row2,col2) || best_clust(row2,col1))
                      {
                        best_clust(row2,col2) = 0;
                        best_clust(row2,col1) = 1;
                      }
                    }
                  }
                }
              }
            }
            unique_cols(best_clust);
            sort_mat(best_clust, col1 + 1);
          }
          best_clust = remove_empty_columns(best_clust);
          best_size = best_clust.n_cols;
        } else {
          for(unsigned int col1 = 0; col1 < best_clust.n_cols; col1++)
          {
            for(unsigned int row = 0; row < best_clust.n_rows; row++)
            {
              if (best_clust(row,col1))
              {
                for(unsigned int col2 = col1 + 1; col2 < best_clust.n_cols; col2++)
                {
                  best_clust(row,col2) = 0;
                }
              }
            }
            unique_cols(best_clust);
            sort_mat(best_clust, col1 + 1);
          }
          best_clust = remove_empty_columns(best_clust);
          best_size = best_clust.n_cols;
        }
      }

      nb_clust += best_size;

      arma::ivec clust(y.n_rows);
      arma::mat cntr(best_size, gr_n_cols + 1);
      arma::mat err(best_size, gr_n_cols);
      std::vector<arma::mat> s(best_size);
      
      for(unsigned int c = 0; c < best_size; c++)
      {
        unsigned int clust_size = arma::accu(best_clust.col(c));
        arma::mat tmp(clust_size, gr_n_cols);
        int j = 0;
        for(unsigned int i = 0; i < best_clust.n_rows; i++)
        {
          if(best_clust(i,c))
          {
            tmp.row(j) = gr.row(i);
            j++;
            clust(i) = k;
          }
        }
        cntr.row(c).cols(0,gr_n_cols - 1) = arma::mean(tmp, 0);
        cntr.row(c).col(gr_n_cols) = clust_size;
        err.row(c) = arma::stddev(tmp, 0);
        s[c] = tmp;
        k++;
      }
      
      clusters[index] = clust;
      centers[index] = cntr;
      errors[index] = err;
      std::vector<int> sizes(best_size);
      for(unsigned int i = 0; i < best_size; i++)
      {
        sizes[i] = arma::accu(best_clust.col(i));
      }
      // arma::rowvec qual(gr_n_cols, arma::fill::zeros);
      // for(unsigned int c = 0; c < best_size; c++)
      // {
      //   if(sizes[c] == 1)
      //   {
      //     continue;
      //   }
      //   for(unsigned int col = 0; col <  s[c].n_cols; col++)
      //   {
      //     double si = 0;
      //     for(unsigned int i = 0; i <  s[c].n_rows; i++)
      //     {
      //       double a = 0;
      //       double b = std::numeric_limits<double>::max_digits10;
      //       for(unsigned int j = 0; j <  s[c].n_rows; j++)
      //       {
      //         a += std::pow(s[c](i, col) - s[c](j, col), 2);
      //       }
      //       a = 1./(sizes[c] - 1.) * std::sqrt(a);
      //       for(unsigned int v = 0; v < best_size; v++)
      //       {
      //         if(v == c)
      //         {
      //           continue;
      //         }
      //         double tmp = 0.;
      //         for(unsigned int j = 0; j <  s[v].n_rows; j++)
      //         {
      //           tmp += std::pow(s[c](i, col) - s[v](j, col), 2);
      //         }
      //         tmp = 1./sizes[v] * std::sqrt(tmp);
      //         b = std::min(b, tmp);
      //       }
      //       if(std::max(a, b) == 0) {
      //         si += 0;
      //       } else {
      //         si += (b - a)/std::max(a, b);
      //       }
      //     }
      //     qual(col) += si;
      //   }
      // }
      // qual = qual/gr.n_rows;
      // quality[index] = qual;
      index++;
    }
    // else  // row-wise
    // {
    //   std::vector<arma::uvec> clusts(y.n_rows);
    //   dist_rowwise(y, clusts, radius[x], method, start, end);
    //   arma::ivec nums(y.n_rows);
    //   for(unsigned int i = 1; i <= y.n_rows; i++) {
    //     nums(i - 1) = i;
    //   }
    //   int max = nums.size() * (nums.size() + 1) / 2;
    //   
    //   std::sort(clusts.begin(), clusts.end(), [&](const arma::uvec & a, const arma::uvec & b) -> bool 
    //   { 
    //     int lhs = max * arma::accu(a) + arma::accu(a % nums);
    //     int rhs = max * arma::accu(b) + arma::accu(b % nums);
    //     return lhs > rhs;
    //   }) ;     
    //   
    //   clusts.erase(std::unique(clusts.begin(), clusts.end(), [](const arma::uvec & a, const arma::uvec & b) -> bool{
    //     return arma::all(a == b);
    //   }), clusts.end());
    //   
    //   
    //   std::vector<arma::uvec> clusterMap;
    //   
    //   makeClusterMap(clusts, clusterMap);
    //   
    //   nb_clust += clusterMap.size();
    //   
    //   
    //   arma::ivec clust(y.n_rows);
    //   arma::mat cntr(clusterMap.size(), groups[x].n_cols);
    //   arma::mat err(clusterMap.size(), groups[x].n_cols);
    //   std::vector<arma::mat> s(clusterMap.size());
    //   for(unsigned int c = 0; c < clusterMap.size(); c++)
    //   {
    //     arma::mat tmp(arma::accu(clusterMap[c]), groups[x].n_cols);
    //     int j = 0;
    //     for(unsigned int i = 0; i < clusterMap[c].size(); i++)
    //     {
    //       if(clusterMap[c](i))
    //       {
    //         tmp.row(j) = groups[x].row(i);
    //         j++;
    //         clust(i) = k;
    //       }
    //     }
    //     cntr.row(c) = arma::mean(tmp, 0);
    //     err.row(c) = arma::stddev(tmp, 0);
    //     s[c] = tmp;
    //     k++;
    //   }
    //   clusters[x] = clust;
    //   centers[x] = cntr;
    //   errors[x] = err;
    //   
    //   std::vector<int> sizes(clusterMap.size());
    //   
    //   for(unsigned int i = 0; i < clusterMap.size(); i++)
    //   {
    //     sizes[i] = arma::accu(clusterMap[i]);
    //   }
    //   
    //   arma::rowvec qual(groups[x].n_cols, arma::fill::zeros);
    //   for(unsigned int c = 0; c < clusterMap.size(); c++)
    //   {
    //     if(sizes[c] == 1)
    //     {
    //       continue;
    //     }
    //     for(unsigned int col = 0; col <  s[c].n_cols; col++)
    //     {
    //       double si = 0;
    //       for(unsigned int i = 0; i <  s[c].n_rows; i++)
    //       {
    //         double a = 0;
    //         double b = std::numeric_limits<double>::max_digits10;
    //         for(unsigned int j = 0; j <  s[c].n_rows; j++)
    //         {
    //           a += std::pow(s[c](i, col) - s[c](j, col), 2);
    //         }
    //         a = 1./(sizes[c] - 1.) * std::sqrt(a);
    //         for(unsigned int v = 0; v < clusterMap.size(); v++)
    //         {
    //           if(v == c)
    //           {
    //             continue;
    //           }
    //           double tmp = 0.;
    //           for(unsigned int j = 0; j <  s[v].n_rows; j++)
    //           {
    //             tmp += std::pow(s[c](i, col) - s[v](j, col), 2);
    //           }
    //           tmp = 1./sizes[v] * std::sqrt(tmp);
    //           b = std::min(b, tmp);
    //         }
    //         si += (b - a)/std::max(a, b);
    //       }
    //       qual(col) += si;
    //     }
    //   }
    //   qual = qual/groups[x].n_rows;
    //   quality[x] = qual;
    // }
    if (last_invoke_ != time(nullptr) && verbose)
    {
      last_invoke_ = time(nullptr);
      Rcpp::Rcout << std::fixed << std::setprecision(1) << rows_done * 100./m.n_rows << "%     \r"; 
    }
  }
  if (verbose)
  {
    Rcpp::Rcout << std::fixed << std::setprecision(1) << 100. << "%     " << std::endl;
  }
  std::vector<int> clusters_vec(m.n_rows);
  int c = 0;
  for(unsigned int i = 0; i < clusters.size(); i++)
  {
    for(unsigned int j = 0; j < clusters[i].size(); j++)
    {
      clusters_vec[c] = clusters[i][j];
      c++;
    }
  }
  
  c = 0;
  arma::mat centers_mat(nb_clust, gr_n_cols + 1);
  arma::mat errors_mat(nb_clust, gr_n_cols);
  for(unsigned int i = 0; i < centers.size(); i++)
  {
    for(unsigned int j = 0; j < centers[i].n_rows; j++)
    {
      centers_mat.row(c) = centers[i].row(j);
      errors_mat.row(c) = errors[i].row(j);
      c++;
    }
  }
  // Rcpp::Rcout << "here";
  // arma::mat quality_mat(quality.size(), gr_n_cols);
  // for(unsigned int i = 0; i < quality.size(); i++)
  // {
  //   if(i == (quality.size() - 1))
  //   quality_mat.row(i) = quality[i];
  // }

  
  
  return Rcpp::List::create(Rcpp::_["clusters"] = clusters_vec, Rcpp::_["centers"] = centers_mat, Rcpp::_["SD"] = errors_mat);//, Rcpp::_["quality"] = quality_mat);
}





