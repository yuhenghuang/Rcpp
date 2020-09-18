#include <Rcpp.h>

// [[Rcpp:plugins(cpp11)]]

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
DateVector test_date() {
  Date d1("2020-01-01");
  Date d2("2020-02-02", "%Y-%m-%d");
  Date d3(3); // "1970-01-01" + days(3);
  Date d4(2020, 1, 2);

  Rcout << d2 - d1 << "\n";
  Rcout << (d2 > d1) << endl;

  DateVector dates(2);
  dates[0] = d1 + 1;
  dates[1] = d4;

  Date d("2016-1-1");
  Rcout << d.getDay() << endl;     //1
  Rcout << d.getMonth() << endl;   //1
  Rcout << d.getYear() << endl;    //2016
  Rcout << d.getWeekday() << endl; //6
  Rcout << d.getYearday() << endl; //1

  return dates;
}


// [[Rcpp::export]]
DatetimeVector test_dttm() {
  // timezone is UTC
  Datetime dt1(10.1); // "1970-01-01 00:00:00 UTC" + 10.1 sec

  // timezone is local timezone
  Datetime dt2("2000-01-01T00", "%Y-%m-%dT%H");

  DatetimeVector dttms(2);
  dttms[0] = dt1;
  dttms[1] = dt2;

  return dttms;
}
