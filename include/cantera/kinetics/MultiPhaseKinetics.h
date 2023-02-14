/**
 * @file MultiPhaseKinetics.h
 * @ingroup chemkinetics
 */

#ifndef MULTIPHASE_KINETICS_H
#define MULTIPHASE_KINETICS_H

#include "Kinetics.h"
#include "InterfaceKinetics.h"

namespace Cantera
{

class MultiPhaseKinetics : public InterfaceKinetics
{
public:
  MultiPhaseKinetics(ThermoPhase* thermo = 0) :
    InterfaceKinetics(thermo)
  {}

  virtual ~MultiPhaseKinetics() {}

  // init is called after adding all phases
  virtual void init() override {
    std::cout << "=== I'm initing ==" << std::endl;
    InterfaceKinetics::init();
    
    // vapor phase is surface phase + 1
    size_t ksurf = reactionPhaseIndex();
    m_vapor_phase = &thermo(ksurf + 1);
  }

  virtual void updateROP() override;

  // equilibrium constants for reverse reactions
  virtual void getEquilibriumConstants(doublereal* kc) {
    for (size_t i = 0; i < nReactions(); i++) {
      kc[i] = 1.;
    }
  }

  virtual void getActivityConcentrations(doublereal* const conc) {
    updateROP();
    copy(m_actConc.begin(), m_actConc.end(), conc);
  }

protected:
  //! pointer to the gas phase
  ThermoPhase const* m_vapor_phase;
};

}

#endif
