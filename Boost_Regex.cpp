#include <Rcpp.h>

#include <string>
#include <boost/regex.hpp>

using namespace Rcpp;

bool validate_card_format(const std::string &s) {
  static const boost::regex e("(\\d{4}[- ]){3}\\d{4}");
  return boost::regex_match(s, e);
}

const boost::regex e("\\A(\\d{3,4})[- ]?(\\d{4})[- ]?(\\d{4})[- ]?(\\d{4})\\z");
const std::string machine_format("\\1\\2\\3\\4");
const std::string human_format("\\1-\\2-\\3-\\4");

std::string machine_readable_card_number(const std::string& s) {
  return boost::regex_replace(s, e, machine_format, boost::match_default | boost::format_sed);
}

std::string human_readable_card_number(const std::string& s) {
  return boost::regex_replace(s, e, human_format, boost::match_default | boost::format_sed);
}

// [[Rcpp::export]]
DataFrame regex_demo(const StringVector &s) {
  int n = s.size();

  LogicalVector valid(n);
  StringVector machine(n), human(n);

  std::string temp;
  for (int i=0; i<n; ++i) {
    temp = s[i];
    valid[i] = validate_card_format(temp);
    machine[i] = machine_readable_card_number(temp);
    human[i] = human_readable_card_number(temp);
  }

  return DataFrame::create(
    _["input"] = s,
    _["valid"] = valid,
    _["machine"] = machine,
    _["human"] = human
  );
}