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
    std::mutex m; // only one thread can either uses push or pop.
    std::condition_variable cv;
    T end_condition;
    Predicate pred;

  public:
    ThreadQueue(const T& ec): end_condition(ec) {}

    void push(T&& elem) {
      std::lock_guard<std::mutex> lock(m);
      
      q.push(elem);

      cv.notify_one();
    }

    bool pop(T& elem) {
      // own mutex
      std::unique_lock<std::mutex> lock(m);
      /*
        1. mutex released and execution suspended before awaking (being awakened)
        2. Being awakened by notify_one() ( or notify_all() )
        3. mutex reacquired, check the condition
        4. If true -> move to next line, else -> go back to 1
      */
      cv.wait(lock, [this](){ return !q.empty(); });

      // if is not necessary?
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
