#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

namespace {

arma::mat as_optional_mat(const Rcpp::Nullable<Rcpp::NumericMatrix>& maybe_mat,
                          const std::size_t n_rows) {
  if (maybe_mat.isNull()) {
    return arma::mat(n_rows, 0, arma::fill::zeros);
  }
  return Rcpp::as<arma::mat>(maybe_mat);
}

int first_bad_row(const arma::mat& pp,
                  const arma::mat& mask,
                  const double tol_xi) {
  if (pp.n_cols == 0 || mask.n_rows == 0) {
    return 0;
  }

  for (arma::uword ii = 0; ii < mask.n_rows; ++ii) {
    arma::uvec observed = arma::find(mask.row(ii) > 0.5);
    if (observed.n_elem == 0) {
      continue;
    }

    arma::mat sub_pp = pp.rows(observed);
    arma::mat gram = sub_pp.t() * sub_pp;
    if (arma::rcond(gram) < tol_xi) {
      return static_cast<int>(ii + 1);
    }
  }

  return 0;
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List pls_component_step_cpp(const Rcpp::NumericMatrix& xxwotna_r,
                                  const Rcpp::NumericMatrix& xxna_r,
                                  const Rcpp::NumericVector& tempww_r,
                                  Rcpp::Nullable<Rcpp::NumericMatrix> prev_pp_r = R_NilValue,
                                  Rcpp::Nullable<Rcpp::NumericMatrix> predict_na_r = R_NilValue,
                                  const double tol_xi = 1e-12,
                                  const bool check_xx = false,
                                  const bool check_predict = false) {
  arma::mat xxwotna = Rcpp::as<arma::mat>(xxwotna_r);
  arma::mat xxna = Rcpp::as<arma::mat>(xxna_r);
  arma::vec tempww = Rcpp::as<arma::vec>(tempww_r);

  double tempww_norm = std::sqrt(arma::dot(tempww, tempww));
  arma::vec tempwwnorm;
  if (tempww_norm > 0) {
    tempwwnorm = tempww / tempww_norm;
  } else {
    tempwwnorm = arma::zeros<arma::vec>(tempww.n_elem);
  }

  arma::vec denom_t = xxna * arma::square(tempwwnorm);
  arma::vec temptt = (xxwotna * tempwwnorm) / denom_t;

  arma::vec denom_p = xxna.t() * arma::square(temptt);
  arma::vec temppp = (xxwotna.t() * temptt) / denom_p;
  arma::mat residxx = xxwotna - temptt * temppp.t();

  arma::mat prev_pp = as_optional_mat(prev_pp_r, xxwotna.n_cols);
  arma::mat combined_pp(xxwotna.n_cols, prev_pp.n_cols + 1, arma::fill::zeros);
  if (prev_pp.n_cols > 0) {
    combined_pp.cols(0, prev_pp.n_cols - 1) = prev_pp;
  }
  combined_pp.col(combined_pp.n_cols - 1) = temppp;

  int bad_xx_row = 0;
  if (check_xx) {
    bad_xx_row = first_bad_row(combined_pp, xxna, tol_xi);
  }

  int bad_predict_row = 0;
  if (check_predict) {
    arma::mat predict_na = as_optional_mat(predict_na_r, 0);
    bad_predict_row = first_bad_row(combined_pp, predict_na, tol_xi);
  }

  return Rcpp::List::create(
    Rcpp::Named("tempwwnorm") = tempwwnorm,
    Rcpp::Named("temptt") = temptt,
    Rcpp::Named("temppp") = temppp,
    Rcpp::Named("residXX") = residxx,
    Rcpp::Named("bad_xx_row") = bad_xx_row,
    Rcpp::Named("bad_predict_row") = bad_predict_row
  );
}
