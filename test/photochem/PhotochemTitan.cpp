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

    // set wavelength
    vector<double> wavelength(10);
    vector<double> actinic_flux(10);

    for (int i = 0; i < 10; i++) {
      wavelength[i] = 20.0 + i * 20.0;
      actinic_flux[i] = 1.0;
    }

    kin->setWavelength(wavelength.data(), wavelength.size());
    kin->updateActinicFlux(actinic_flux.data());
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
  ASSERT_EQ(kin->nWavelengths(), 10);
}

TEST_F(PhotochemTitan, check_fwd_rate_constants) {
  vector<double> kfwd(kin->nReactions());

  kin->getFwdRateConstants(kfwd.data());

  ASSERT_NEAR(kfwd[0], 3.06820e-14, 1.0e-18);
  ASSERT_NEAR(kfwd[1], 3.2e-16, 1.0e-18);
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
