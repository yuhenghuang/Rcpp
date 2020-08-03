#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List initialize_list(StringVector &names,
                     IntegerVector &ids) {
  // ...
  unsigned int n = names.size();
  List res(n);
  Environment env_g = Environment::global_env();
  Environment env_p = env_g["Person"];
  Function new_person = env_p["new"];

  for (unsigned int i=0; i<n; ++i)
    res[i] = new_person(String(names[i]), int(ids[i]));
  return res;
}

// [[Rcpp::export]]
String access_method(std::string foo_name) {
  Environment g_env = Environment::global_env();
  Environment p_env = g_env["Person"];
  Function new_person = p_env["new"];
  Environment new_p;
  String name = "Jake";
  int id = 1;
  
  new_p = new_person(name, id);
  Function foo = new_p[foo_name];
  String res = foo();
  
  return res;
}

// [[Rcpp::export]]
List initialize_list_pkg(unsigned int size){
  /*
  For classes in packages.
  Most common case.
  */
  List res(size);
  Environment package_env("package:some_package");
  Environment class_env = package_env["some_class"];
  Function new_instance = class_env["new"];
  
  for(unsigned int i = 0; i < size; i++){
    Environment new_i;
    new_i = new_instance();
    res[i] = new_i;
  }
  
  return res;
}

// [[Rcpp::export]]
void give_item(Environment &person,
               String &item) {
  // ...
  Function give_i = person["give_item"];
  give_i(item);
}