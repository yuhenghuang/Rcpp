#include <Rcpp.h>

using namespace Rcpp;


class Uniform {
  private:
    double min, max;
  public:
    Uniform(double _min, double _max): 
      min(_min), max(_max) { }

    NumericVector draw(int n) {
      if (n<1) return NumericVector();
      RNGScope scope;
      return runif(n, min, max);
    }

    const double &getmin() const {
      return min;
    }

    const double &getmax() const {
      return max;
    }
};




RCPP_MODULE(unif_module) {

  Rcpp::class_<Uniform>("Uniform")

  .constructor<double, double>()

  .method("draw", &Uniform::draw, "Draw n random values following the distribution")


  .property("min", &Uniform::getmin)
  .property("max", &Uniform::getmax)
  ;

}


class Example {
  public:
    Example(int x_, int y_): x(x_), y(y_) {
      Rcout << __PRETTY_FUNCTION__ << "\n";
    }

    Example(std::string x_, std::string y_): x(x_.size()), y(y_.size()) {
      Rcout << __PRETTY_FUNCTION__ << "\n";
    }

    int add() const { 
      return x + y; 
    }

private:
    int x, y;
};


template <typename T>
bool example_validator(SEXP* args, int nargs) {
  const int RTYPE = Rcpp::traits::r_sexptype_traits<T>::rtype;
  
  for (int i=0; i<nargs; ++i)
    if (TYPEOF(args[i])!=RTYPE)
      return false;
  return true;
}

/*
template <typename T>
bool universal_validator(SEXP x) {
  const int RTYPE = Rcpp::traits::r_sexptype_traits<T>::rtype;
  return TYPEOF(x) == RTYPE;
}
*/

template <typename... Types>
bool universal_validator(SEXP* args, int nargs) {
  return universal_validator<Types...>(args, nargs, 0);
}

template <typename T = void, typename... Types>
bool universal_validator(SEXP* args, int nargs, int idx) {
  return Rcpp::is<T>(args[idx]) && universal_validator<Types...>(args, nargs, idx+1);
}

template <>
bool universal_validator<>(SEXP* args, int nargs, int idx) {
  return nargs == idx;
}



RCPP_MODULE(example_module) {

  Rcpp::class_<Example>("Example")
  
  // .constructor<int, int>("(int, int) constructor", example_validator<int>)
  // .constructor<std::string, std::string> ("(string, string) constructor", example_validator<std::string>)
  .constructor<int, int>("(int, int) constructor", universal_validator<int, int>)
  .constructor<std::string, std::string> ("(string, string) constructor", universal_validator<std::string, std::string>)
  .method("add", &Example::add)
  ;
}


