#ifndef PHOTOLYSIS_RATE_HPP
#define PHOTOLYSIS_RATE_HPP

#include "cantera/kinetics/ReactionData.h"
#include "ReactionRate.h"
#include "MultiRate.h"

namespace Cantera
{

struct PhotolysisData: public ReactionData
{
  virtual bool update(const ThermoPhase& phase, const Kinetics& kin) override;
  using ReactionData::update;
};

class PhotolysisRate : public ReactionRate
{
public:
  PhotolysisRate(const AnyMap& node, const UnitStack& rate_units = {})
  {
    setParameters(node, rate_units);
  }

  unique_ptr<MultiRateBase> newMultiRate() const {
    return unique_ptr<MultiRateBase>(new MultiRate<PhotolysisRate, PhotolysisData>);
  }

  const std::string type() const { return "photolysis"; }

  void setParameters(const AnyMap& node, const UnitStack& rate_units);

  void getParameters(AnyMap& node) const;

  void updateFromStruct(const FreezingData& shared_data);

  double evalFromStruct(const FreezingData& shared_data);

  double ddTScaledFromStruct(const FreezingData& shared_data) const;
};

}

#endif
