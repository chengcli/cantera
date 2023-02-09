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
    ready(false), pressure(-1)
  {}

  void update(double T) override {
    throw CanteraError("CondensationData::update",
        "Missing state information: 'CondensationData' requires pressure.");
  }

  void update(double T, double P) override {
    ReactionData::update(T);
    pressure = P;
  }

  using ReactionData::update;

  bool update(const ThermoPhase& phase, const Kinetics& kin) override;

  void resize(size_t nSpecies, size_t nReactions, size_t nPhases) override {
    mole_fractions.resize(nSpecies, 0.);
    ready = true;
  }

  void restore() override {
    ReactionData::restore();

    // only restore if there is a valid buffered value
    if (m_pressure_buf < 0.) {
        return;
    }
    update(temperature, m_pressure_buf);
    m_pressure_buf = -1.;
  }

  virtual void invalidateCache() override {
    ReactionData::invalidateCache();
    pressure = NAN;
  }

  void perturbPressure(double deltaP)
  {
    if (m_pressure_buf > 0.) {
      throw CanteraError("CondensationData::perturbPressure",
      "Cannot apply another perturbation as state is already perturbed.");
    }
    m_pressure_buf = pressure;
    update(temperature, pressure * (1. + deltaP));
  }

  bool      ready; //!< boolean indicating whether vectors are accessible
  double    pressure;
  vector_fp mole_fractions;

protected:
  double m_pressure_buf; //!< buffered pressure
};

class CondensationRate : public ReactionRate
{
public:
  CondensationRate(const AnyMap& node, const UnitStack& rate_units = {}):
    m_vapor_index(-1), m_A(0.), m_p(0.), m_a(0.), m_b(0.)
  {
    setParameters(node, rate_units);
  }

  unique_ptr<MultiRateBase> newMultiRate() const override {
    return unique_ptr<MultiRateBase>(new MultiRate<CondensationRate, CondensationData>);
  }

  const std::string type() const override { return "condensation"; }

  void setParameters(const AnyMap& node, const UnitStack& rate_units) override;

  void setContext(const Reaction& rxn, const Kinetics& kin) override;

  double evalFromStruct(const CondensationData& shared_data);

  // general saturation vapor pressure expressed as
  double saturation_vapor_pressure(double T) {
    return m_p*exp(m_a - m_b/T);
  }

protected:
  void getParameters(AnyMap& node) const override;

  int    m_vapor_index;
  double m_A;  // condensation rate
  double m_p;
  double m_a;
  double m_b;
};

}

#endif
