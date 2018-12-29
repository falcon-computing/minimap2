#ifndef FALCONLIC_FLEXLICENSE_H
#define FALCONLIC_FLEXLICENSE_H
#include <string>
#include <unordered_map>

#include "flexlm/lm_attr.h"
#include "flexlm/lmclient.h"

#include "license.h"

namespace falconlic {
namespace flexlm {

class License {
 public:
  License();
  ~License();
  void add_feature(Feature f);
  int verify();

 private:
  int checkout(Feature f);
  int checkin(Feature f);

  int flex_init();
  void flex_clean();

  std::unordered_map<int, bool> feat_table_;

  VENDORCODE code;
  LM_HANDLE *lm_job = 0;
  struct flexinit_property_handle *initHandle = 0;
};

} // namespace flexlm
} // namespace falconlic
#endif
