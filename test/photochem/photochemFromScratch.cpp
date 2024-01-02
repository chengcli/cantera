#include "gtest/gtest.h"

using namespace Cantera;

class PhotochemFromScratch: public testing::Test {
 public:
  PhotochemFromScratch() {
    // set the initial state
    p1.setState_TPX(1200.0, 5.0*OneAtm, "H2:1.0, O2:1.0");
    p1.equilibrate("TP");
    p1.save("h2o2.xml","tran",false);
    p1.save("h2o2.cti","cti",false);
  }

  shared_ptr<ThermoPhase> p1;
  BulkKinetics kin;
};

TEST_F(PhotochemFromScratch, build) {
  IdealGasMix p2("h2o2.xml","tran");
  p2.setState_TPX(1200.0, 5.0*OneAtm, "H2:1.0, O2:1.0");
  p2.equilibrate("TP");
  p2.save("h2o2.xml","tran",false);
  p2.save("h2o2.cti","cti",false);
}
