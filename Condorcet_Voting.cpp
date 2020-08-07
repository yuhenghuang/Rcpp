#include <RcppArmadillo.h>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using arma::uword;

arma::Mat<int> pair_count_cpp(const arma::Mat<int> &mat_ballot) {
  // ...
  uword num_candidates = mat_ballot.n_rows, num_ballots = mat_ballot.n_cols;
  arma::Mat<int> pairwise(num_candidates, num_candidates);

  int count;
  for (uword i=0; i<num_candidates; ++i) {
    pairwise(i, i) = 0;
    for (uword j=0; j<i; ++j) {
      count = 0;
      for (uword k=0; k<num_ballots; ++k)
        count += mat_ballot(i, k) < mat_ballot(j, k) ? 1 : 0;
      pairwise(i, j) = count;
      pairwise(j, i) = int(num_ballots) - count;
    }
  }

  return pairwise;
}

arma::Mat<int> schulze_cpp(const arma::Mat<int> &pairwise) {
  uword n = pairwise.n_rows;
  arma::Mat<int> schulze(n, n);

  bool flag;
  for (uword i=0; i<n; ++i) {
    schulze(i, i) = 0;
    for (uword j=0; j<i; ++j) {
      flag = pairwise(i, j) > pairwise(j, i);
      schulze(i, j) = flag ? pairwise(i, j) : 0;
      schulze(j, i) = flag ? pairwise(j, i) : 0;
    }
  }

  for (uword i=0; i<n; ++i)
    for (uword j=0; j<n; ++j)
      if (i!=j) 
        for (uword k=0; k<n; ++k)
          if (i!=k && j!=k)
            schulze(j, k) = std::max(schulze(j, k), std::min(schulze(j, i), schulze(i, k)));

  return schulze;
}

std::vector<uword> find_winners(const uword &current_size,
                                const arma::Mat<int> &pairwise) {
  // ...
  int count, maximum=0;
  std::vector<uword> winners;
  for (uword i=0; i<current_size; ++i) {
    count = 0;
    for (uword j=0; j<current_size; ++j)
      count += pairwise(i, j) > pairwise(j, i) ? 1 : 0;
    if (count>maximum) {
      winners.clear();
      maximum = count;
    }
    if (count==maximum) 
      winners.push_back(i);
  }
  return winners;
}

// [[Rcpp::export]]
List condorcet_rank(const StringVector &name_candidates, 
                    const arma::Mat<int> &mat_ballot) {
  // ...
  uword num_candidates = mat_ballot.n_rows;

  uword current_rank = 1, current_size = num_candidates;
  StringVector current_names = name_candidates;
  arma::Mat<int> current_ballot = mat_ballot;

  StringVector name, method;
  IntegerVector rank;
  arma::Mat<int> pairwise;

  std::vector<uword> winners;

  uword winner;
  while (current_rank<=num_candidates) {
    pairwise = pair_count_cpp(current_ballot);
    winners = find_winners(current_size, pairwise);

    if (winners.size()==1) {
      winner = winners.back();
      winners.pop_back();

      rank.push_back(current_rank);
      name.push_back(current_names[winner]);
      method.push_back("Condorcet");

      ++current_rank;
      --current_size;
      current_ballot.shed_row(winner);
      current_names.erase(winner);
    }
    else {
      pairwise = schulze_cpp(pairwise);
      winners = find_winners(current_size, pairwise);
      
      int temp_size = winners.size();
      // must loop backward over winners
      while (!winners.empty()) {
        winner = winners.back();
        winners.pop_back();

        rank.push_back(current_rank);
        name.push_back(current_names[winner]);
        method.push_back("Schulze");

        --current_size;
        current_ballot.shed_row(winner);
        current_names.erase(winner);
      }
      current_rank += temp_size;
    }
  }

  return List::create(
    _["Name"] = name,
    _["Rank"] = rank,
    _["Method"] = method
  );
}