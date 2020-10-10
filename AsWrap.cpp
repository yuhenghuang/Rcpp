#include <Rcpp.h>
#include <boost/numeric/ublas/vector.hpp>


namespace Rcpp {
  
  namespace traits {

    template <typename T> 
    SEXP wrap(const boost::numeric::ublas::vector<T>& obj) {
      const int RTYPE = Rcpp::traits::r_sexptype_traits<T>::rtype;
      // Rcpp::traits::storage_type<RTYPE>::type

      return Rcpp::Vector<RTYPE>(obj.begin(), obj.end());
    }

    template <typename T>
    class Exporter<boost::numeric::ublas::vector<T>> {
      typedef typename boost::numeric::ublas::vector<T> out_type;

      const static int RTYPE = Rcpp::traits::r_sexptype_traits<T>::rtype;
      Rcpp::Vector<RTYPE> vec;

      public:
        Exporter(SEXP x): vec(x) {
          if (TYPEOF(x) != RTYPE)
            throw std::invalid_argument("Wrong R type");
        }

        out_type get() {
          out_type x(vec.size());
          std::copy(vec.begin(), vec.end(), x.begin());
          return x;
        }
    };

  }

}