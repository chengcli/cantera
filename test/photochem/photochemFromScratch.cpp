#include "gtest/gtest.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/kinetics/KineticsFactory.h"

namespace Cantera
{

class PhotochemTitan: public testing::Test {
 public:
  // data
  shared_ptr<ThermoPhase> phase;
  shared_ptr<Kinetics> kin;

  // constructor
  PhotochemTitan() {
    phase = newThermo("../data/ch4_photolysis.yaml");
    kin = newKinetics({phase}, "../data/ch4_photolysis.yaml");

    // set the initial state
    string X = "CH4:0.02 N2:0.98";
    phase->setState_TPX(200.0, OneAtm, X);
  }
};

TEST_F(PhotochemTitan, check_phase) {
  ASSERT_EQ(phase->nElements(), 3);
  ASSERT_EQ(phase->nSpecies(), 8);
}

TEST_F(PhotochemTitan, check_kinetics) {
  ASSERT_EQ(kin->nReactions(), 2);
  ASSERT_EQ(kin->nTotalSpecies(), 8);
  ASSERT_EQ(kin->nPhases(), 1);
}

} // namespace Cantera

int main(int argc, char** argv)
{
  printf("Running main() from PhotochemTitan.cpp\n");
  Cantera::make_deprecation_warnings_fatal();
  Cantera::printStackTraceOnSegfault();
  testing::InitGoogleTest(&argc, argv);
  int result = RUN_ALL_TESTS();
  Cantera::appdelete();
  return result;
}
