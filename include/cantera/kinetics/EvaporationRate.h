#ifndef EVAPORATION_RATE_HPP
#define EVAPORATION_RATE_HPP

#include "cantera/kinetics/ReactionData.h"
#include "ReactionRate.h"
#include "MultiRate.h"

namespace Cantera
{

struct EvaporationData: public ReactionData
{
  EvaporationData()
  {}

  virtual bool update(const ThermoPhase& phase, const Kinetics& kin) override;
  using ReactionData::update;
};

class EvaporationRate : public ReactionRate
{
public:
  EvaporationRate(const AnyMap& node, const UnitStack& rate_units = {}):
    m_vapor_index(-1), m_solid_index(-1), m_A(0.), m_P3(0.), m_T3(1.), m_beta(0.), m_gamma(0.), m_ice(false)
  {
    setParameters(node, rate_units);
  }

  unique_ptr<MultiRateBase> newMultiRate() const override {
    return unique_ptr<MultiRateBase>(new MultiRate<EvaporationRate, EvaporationData>);
  }

  const std::string type() const { return "evaporation"; }

  void setParameters(const AnyMap& node, const UnitStack& rate_units);

  void setContext(const Reaction& rxn, const Kinetics& kin) override;

  void getActivityConcentration(double *actConc, double const* conc,
    EvaporationData const& shared_data);

  double evalFromStruct(const EvaporationData& shared_data);

  // ideal saturation vapor pressure
  double saturation_vapor_pressure(double T) {
    return m_P3*exp((1. - m_T3/T)*m_beta - m_gamma*log(T/m_T3));
  }

protected:
  void getParameters(AnyMap& node) const;

  int    m_vapor_index;
  int    m_solid_index;
  double m_A;  // evaporation rate
  double m_P3;
  double m_T3;
  double m_beta;
  double m_gamma;
  bool   m_ice;
};

}

#endif
