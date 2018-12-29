#ifndef FALCON_LIC_CASE_H
#define FALCON_LIC_CASE_H

#include "license.h"

// built-in license check for genomics products
inline int license_verify() {
#ifdef USELICENSE
  namespace fc   = falconlic;
  namespace fclm = falconlic::flexlm;
#ifdef DEPLOY_aws
  fc::enable_aws();
#endif
#ifdef DEPLOY_hwc
  fc::enable_hwc();
#endif

#ifndef NDEBUG
  fc::set_verbose(3);
#endif

  fclm::enable();
  fclm::add(fclm::FALCON_DNA);
  return fc::license_verify();
#else
  return 0;
#endif
}

#endif
