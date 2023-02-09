#ifndef AUTOCONVERSION_RATE_HPP
#define AUTOCONVERSION_RATE_HPP

#include "cantera/kinetics/ReactionData.h"
#include "ReactionRate.h"
#include "MultiRate.h"

namespace Cantera
{

struct AutoconversionData: public ReactionData
{
  AutoconversionData()
  {}

  bool update(const ThermoPhase& phase, const Kinetics& kin) override;
  using ReactionData::update;
};

class AutoconversionRate : public ReactionRate
{
public:
  // access data
  double autoconversion_rate;

  AutoconversionRate(const AnyMap& node, const UnitStack& rate_units = {})
  {
    setParameters(node, rate_units);
  }

  unique_ptr<MultiRateBase> newMultiRate() const {
    return unique_ptr<MultiRateBase>(new MultiRate<AutoconversionRate, AutoconversionData>);
  }

  const std::string type() const { return "autoconversion"; }

  void setParameters(const AnyMap& node, const UnitStack& rate_units);

  double evalFromStruct(const AutoconversionData& shared_data);

protected:
  void getParameters(AnyMap& node) const;
};

}

#endif
