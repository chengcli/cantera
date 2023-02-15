#ifndef COAGULATION_RATE_HPP
#define COAGULATION_RATE_HPP

#include "cantera/kinetics/ReactionData.h"
#include "ReactionRate.h"
#include "MultiRate.h"

namespace Cantera
{

struct CoagulationData: public ReactionData
{
  CoagulationData()
  {}

  bool update(const ThermoPhase& phase, const Kinetics& kin) override;
  using ReactionData::update;
};

class CoagulationRate : public ReactionRate
{
public:
  // access data
  double autoconversion_rate;

  CoagulationRate(const AnyMap& node, const UnitStack& rate_units = {})
  {
    setParameters(node, rate_units);
  }

  unique_ptr<MultiRateBase> newMultiRate() const {
    return unique_ptr<MultiRateBase>(new MultiRate<CoagulationRate, CoagulationData>);
  }

  const std::string type() const { return "autoconversion"; }

  void setParameters(const AnyMap& node, const UnitStack& rate_units);

  double evalFromStruct(const CoagulationData& shared_data);

protected:
  void getParameters(AnyMap& node) const;
};

}

#endif
