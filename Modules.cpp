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
    Example(): x(0), y(0) {
      Rcout << __PRETTY_FUNCTION__ << "\n";
    }

    Example(int x_): x(x_), y(0) {
      Rcout << __PRETTY_FUNCTION__ << "\n";
    }

    Example(int x_, int y_): x(x_), y(y_) {
      Rcout << __PRETTY_FUNCTION__ << "\n";
    }

    Example(std::string x_, std::string y_): x(x_.size()), y(y_.size()) {
      Rcout << __PRETTY_FUNCTION__ << "\n";
    }

    Example(std::string x_, int y_): x(x_.size()), y(y_){
      Rcout << __PRETTY_FUNCTION__ << "\n";
    }

    Example(int x_, int y_, int z_): x(x_), y(y_+z_) {
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
  // return universal_validator<Types...>(args, nargs, 0);
  return universal_validator<Types...>(args, args+nargs);
}

template <typename T = void, typename... Types>
bool universal_validator(SEXP* args, int nargs, int idx) {
  if (idx>=nargs) return false;
  // optional type traits
  typedef typename Rcpp::traits::remove_const_and_reference<T>::type _Tp;
  return Rcpp::is<_Tp>(args[idx]) && universal_validator<Types...>(args, nargs, idx+1);
}

template <>
bool universal_validator<>(SEXP* args, int nargs, int idx) {
  return nargs == idx;
}


template <typename T = void, typename... Types>
bool universal_validator(SEXP* args, SEXP* end) {
  if (args==end) return false;
  typedef typename Rcpp::traits::remove_const_and_reference<T>::type _Tp;
  return Rcpp::is<_Tp>(*args) && universal_validator<Types...>(args+1, end);
}


template <>
bool universal_validator(SEXP* args, SEXP* end) {
  return args == end;
}

RCPP_EXPOSED_CLASS(Example);

#define CTOR(...) \
  constructor<__VA_ARGS__>("(" #__VA_ARGS__ ") constructor", universal_validator<__VA_ARGS__>)

RCPP_MODULE(example_module) {

  Rcpp::class_<Example>("Example")
  
  .constructor()

  // .constructor<int, int>("(int, int) constructor", example_validator<int>)
  // .constructor<std::string, std::string> ("(string, string) constructor", example_validator<std::string>)

  /*
  .constructor<int>("(int) constructor", universal_validator<int>)
  .constructor<int, int>("(int, int) constructor", universal_validator<int, int>)
  .constructor<std::string, std::string> ("(string, string) constructor", universal_validator<std::string, std::string>)
  .constructor<std::string, int> ("(string, int) constructor", universal_validator<std::string, int>)
  .constructor<int, int, int>("(int, int, int) constructor", universal_validator<int, int, int>)
  */

  .CTOR(int)
  .CTOR(int, int)
  .CTOR(std::string, std::string)
  .CTOR(std::string, int)
  .CTOR(int, int, int)
  .method("add", &Example::add)
  ;
}


/*** R

# default constructor
example_def <- new(Example)
example_def$add()

# int
example_int <- new (Example, 2L)
example_int$add()

# int, int
example_int_int <- new(Example, 1L, 3L) 
example_int_int$add()

# str, str
example_str_str <- new(Example, "abc", "de")
example_str_str$add()

# str, int
example_str_int <- new(Example, "abc", 4L)
example_str_int$add()

# int, int, int
example_iii <- new(Example, 1L, 4L, 3L)
example_iii$add()

*/