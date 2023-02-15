/**
 * @file EvaporationKinetics.h
 * @ingroup chemkinetics
 */

#ifndef EVAPORATION_KINETICS_H
#define EVAPORATION_KINETICS_H

#include "Kinetics.h"
#include "InterfaceKinetics.h"

namespace Cantera
{

class EvaporationKinetics : public InterfaceKinetics
{
public:
  EvaporationKinetics(ThermoPhase* thermo = 0) :
    InterfaceKinetics(thermo)
  {}

  virtual ~EvaporationKinetics() {}

  // init is called after adding all phases
  virtual void init() override {
    std::cout << "=== I'm initing evaporation ==" << std::endl;
    InterfaceKinetics::init();
    
    // vapor phase is surface phase + 1
    size_t ksurf = reactionPhaseIndex();
    m_vapor_phase = &thermo(ksurf + 1);
  }

  virtual void updateROP() override;

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
