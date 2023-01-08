#ifndef FREEZING_RATE_HPP
#define FREEZING_RATE_HPP

#include "cantera/kinetics/ReactionData.h"
#include "ReactionRate.h"
#include "MultiRate.h"

namespace Cantera
{

struct FreezingData: public ReactionData
{
  virtual bool update(const ThermoPhase& phase, const Kinetics& kin) override;
  using ReactionData::update;
};

class FreezingRate final : public ReactionRate
{
public:
  // access data
  double freezing_rate;
  double freezing_temperature;

  // functions
  FreezingRate(): freezing_rate(0.), freezing_temperature(0.)
  {}

  FreezingRate(const AnyMap& node, const UnitStack& rate_units={})
  {
      setParameters(node, rate_units);
  }

  unique_ptr<MultiRateBase> newMultiRate() const {
    return unique_ptr<MultiRateBase>(new MultiRate<FreezingRate, FreezingData>);
  }

  //! Identifier of reaction rate type
  const std::string type() const { return "frezing"; }
  
  void setParameters(const AnyMap& node, const UnitStack& rate_units);

  void getParameters(AnyMap& node) const;

  //! Update information specific to reaction
  /*!
   *  @param shared_data  data shared by all reactions of a given type
   */
  void updateFromStruct(const FreezingData& shared_data);

  //! Evaluate reaction rate
  /*!
   *  @param shared_data  data shared by all reactions of a given type
   */
  double evalFromStruct(const FreezingData& shared_data);

  //! Evaluate derivative of reaction rate with respect to temperature
  //! divided by reaction rate
  /*!
   *  @param shared_data  data shared by all reactions of a given type
   */
  double ddTScaledFromStruct(const FreezingData& shared_data) const;

  //! Check basic syntax and settings of reaction rate expression
  void check(const std::string& equation);

  //! Check to make sure that the rate expression is finite over a range of
  //! temperatures at each interpolation pressure. This is potentially an
  //! issue when one of the Arrhenius expressions at a particular pressure
  //! has a negative pre-exponential factor.
  void validate(const std::string& equation, const Kinetics& kin);
};

}

#endif
