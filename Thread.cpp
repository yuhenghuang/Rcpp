#include <Rcpp.h>
#include <unistd.h>
#include <queue>
#include <mutex>
#include <condition_variable>

// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

namespace std {
  template<>
  struct equal_to<string> {
    bool operator()(const string& lhs, const string& rhs) const {
      return lhs.compare(rhs) == 0;
    }
  };
}

template<typename T, class Predicate=std::equal_to<T>>
class ThreadQueue {
  private:
    std::queue<T> q;
    std::mutex m;
    std::condition_variable cv;
    T end_condition;
    Predicate pred;

  public:
    ThreadQueue(const T& ec): end_condition(ec) {}

    void push(T&& elem) {
      std::unique_lock<std::mutex> lck(m);
      q.push(elem);
      cv.notify_one();
    }

    bool pop(T& elem) {
      std::unique_lock<std::mutex> lck(m);
      // why the overloaded version, bool !q.empty(), does not work...
      cv.wait(lck, [this](){ return !q.empty(); });

      if (!q.empty()) {
        elem = std::move(q.front());
        q.pop();
      }

      return pred(elem, end_condition) ? false : true;
    }

};



void read_stdin(ThreadQueue<std::string>& q, std::string end_condition) {
  std::string line;
  while (getline(std::cin, line)) {
    q.push(std::move(line));
  }
  q.push(std::move(end_condition));
}


// [[Rcpp::export]]
void test_stdin() {
  std::string line;
  while (getline(std::cin, line)) {
    sleep(10);
    Rcout << "Out: " << line << std::endl;
  }
}
