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

  int iCH4 = kin->kineticsSpeciesIndex("CH4");
  int iCH3 = kin->kineticsSpeciesIndex("CH3");
  int i1CH2 = kin->kineticsSpeciesIndex("(1)CH2");
  int i3CH2 = kin->kineticsSpeciesIndex("(3)CH2");
  int iCH = kin->kineticsSpeciesIndex("CH");
  int iH2 = kin->kineticsSpeciesIndex("H2");
  int iH = kin->kineticsSpeciesIndex("H");

  ASSERT_EQ(iCH4, 0);
  ASSERT_EQ(iCH3, 1);
  ASSERT_EQ(i1CH2, 2);
  ASSERT_EQ(i3CH2, 3);
  ASSERT_EQ(iCH, 4);
  ASSERT_EQ(iH2, 5);
  ASSERT_EQ(iH, 6);

  double kCH4 = kin->productStoichCoeff(iCH4, 0);
  ASSERT_NEAR(kCH4, 0.635657, 1.0e-4);

  double kCH3 = kin->productStoichCoeff(iCH3, 0);
  ASSERT_NEAR(kCH3, 0.142168, 1.0e-4);

  double k1CH2 = kin->productStoichCoeff(i1CH2, 0);
  ASSERT_NEAR(k1CH2, 0.0978033, 1.0e-4);

  double k3CH2 = kin->productStoichCoeff(i3CH2, 0);
  ASSERT_NEAR(k3CH2, 0.0377844, 1.0e-4);

  double kCH = kin->productStoichCoeff(iCH, 0);
  ASSERT_NEAR(kCH, 0.0865869, 1.0e-4);

  double kH2 = kin->productStoichCoeff(iH2, 0);
  ASSERT_NEAR(kH2, 0.18439, 1.0e-4);

  double kH = kin->productStoichCoeff(iH, 0);
  ASSERT_NEAR(kH, 0.304324, 1.0e-4);

  ASSERT_NEAR(kCH4 + kCH3 + k1CH2 + k3CH2 + kCH, 1.0, 1.0e-14);
  ASSERT_NEAR(4 * kCH4 + 3 * kCH3 + 2 * k1CH2 + 2 * k3CH2 + kCH + 2 * kH2 + kH, 4.0, 1.0e-14);
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
