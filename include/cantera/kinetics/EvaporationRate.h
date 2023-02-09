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
  // access data
  double evaportion_rate;

  EvaporationRate(const AnyMap& node, const UnitStack& rate_units = {})
  {
    setParameters(node, rate_units);
  }

  unique_ptr<MultiRateBase> newMultiRate() const override {
    return unique_ptr<MultiRateBase>(new MultiRate<EvaporationRate, EvaporationData>);
  }

  const std::string type() const { return "evaporation"; }

  void setParameters(const AnyMap& node, const UnitStack& rate_units);

  // void setContext(const Reaction& rxn, const Kinetics& kin) override;

  double evalFromStruct(const EvaporationData& shared_data);

protected:
  void getParameters(AnyMap& node) const;
};

}

#endif
