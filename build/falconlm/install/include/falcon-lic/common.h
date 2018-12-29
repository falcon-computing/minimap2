#ifndef FCS_LICENSE_COMMON_H
#define FCS_LICENSE_COMMON_H

namespace falconlic {
extern int verbose;

void print_error(const char* msg, ...);
void print_warn(const char* msg, ...);
void print_info(const char* msg, ...);

} // namespace falconlic
#endif
