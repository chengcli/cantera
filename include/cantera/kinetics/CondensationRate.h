#ifndef CONDENSATION_RATE_HPP
#define CONDENSATION_RATE_HPP

#include "cantera/kinetics/ReactionData.h"
#include "ReactionRate.h"
#include "MultiRate.h"

namespace Cantera
{

struct CondensationData: public ReactionData
{
  CondensationData():
    ready(false)
  {}

  using ReactionData::update;

  bool update(const ThermoPhase& phase, const Kinetics& kin) override;

  void resize(size_t nSpecies, size_t nReactions, size_t nPhases) override {
    phase_temperatures.resize(nPhases, 0.);
    ready = true;
  }

  bool      ready; //!< boolean indicating whether vectors are accessible
  vector_fp phase_temperatures;
};

class CondensationRate : public ReactionRate
{
public:
  CondensationRate(const AnyMap& node, const UnitStack& rate_units = {}):
    m_vapor_index(-1), m_solid_index(-1), m_A(0.), m_P3(0.), m_T3(1.), m_beta(0.), m_gamma(0.), m_ice(false)
  {
    setParameters(node, rate_units);
  }

  unique_ptr<MultiRateBase> newMultiRate() const override {
    return unique_ptr<MultiRateBase>(new MultiRate<CondensationRate, CondensationData>);
  }

  const std::string type() const override { return "condensation"; }

  void setParameters(const AnyMap& node, const UnitStack& rate_units) override;

  void setContext(const Reaction& rxn, const Kinetics& kin) override;

  void getActivityConcentration(double *actConc, double const* conc,
    CondensationData const& shared_data);

  double evalFromStruct(const CondensationData& shared_data);

  // ideal saturation vapor pressure
  double saturation_vapor_pressure(double T) {
    return m_P3*exp((1. - m_T3/T)*m_beta - m_gamma*log(T/m_T3));
  }

protected:
  void getParameters(AnyMap& node) const override;

  int    m_vapor_index;
  int    m_solid_index;
  double m_A;  // condensation rate
  double m_P3;
  double m_T3;
  double m_beta;
  double m_gamma;
  bool   m_ice;
};

}

#endif
