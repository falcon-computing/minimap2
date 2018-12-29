#ifndef FALCONLIC_LICENSE_H
#define FALCONLIC_LICENSE_H

namespace falconlic {

const int SUCCESS        = 0;
const int API_ERROR      = -1;
const int SIGN_ERROR     = -2;
const int PARSE_ERROR    = -3;
const int PRD_VERF_ERROR = -4;
const int FLEXLM_ERROR   = -5;

namespace flexlm {
  typedef enum {
    FALCON_XILINX = 0, 
    FALCON_ALTERA = 1, 
    FALCON_CUSTOM = 2, 
    FALCON_RT = 3, 
    FALCON_DNA = 4
  } Feature;

  void enable();
  void add(Feature f);
  void clean();
}

void enable_aws();       // default: disabled
void enable_hwc();       // default: disabled

void set_verbose(int v); // default: 0

int license_verify();
int license_clean();

} // namespace falconlic

#endif
