#ifndef FCS_LICENSE_TIMER_H
#define FCS_LICENSE_TIMER_H

#include <sys/time.h>
#include <time.h>

namespace falconlic {
class Timer {
 public: 
  Timer(std::string func = ""): func_(func) {
    start_ts_ = getUs(); 
  }
  ~Timer() {
    print_info("%s takes %ld us\n", func_.c_str(), getUs()-start_ts_);
  }
  
 private:
  inline uint64_t getUs() {
    struct timespec tr;
    clock_gettime(CLOCK_REALTIME, &tr);
    return (uint64_t)tr.tv_sec*1e6 + tr.tv_nsec/1e3;
  }

  std::string func_;
  uint64_t start_ts_;
};

} // namespace falconlic
#endif
