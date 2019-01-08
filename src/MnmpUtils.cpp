#include "MnmpUtils.h"

#include <sys/time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/syscall.h>

/*
double cputime()
{
        struct rusage r;
        getrusage(RUSAGE_SELF, &r);
        return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

double realtime()
{
        struct timeval tp;
        struct timezone tzp;
        gettimeofday(&tp, &tzp);
        return tp.tv_sec + tp.tv_usec * 1e-6;
}
*/

uint32_t getTid() {
  static bool lacks_gettid = false;

  if (!lacks_gettid) {
    pid_t tid = syscall(__NR_gettid);
    if (tid != -1) {
      return (uint32_t)tid;
    }
    // Technically, this variable has to be volatile, but there is a small
    // performance penalty in accessing volatile variables and there should
    // not be any serious adverse effect if a thread does not immediately see
    // the value change to "true".
    lacks_gettid = true;
  }
}
