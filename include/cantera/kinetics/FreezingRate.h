#ifndef FREEZING_RATE_HPP
#define FREEZING_RATE_HPP

#include "cantera/kinetics/ReactionData.h"
#include "ReactionRate.h"
#include "MultiRate.h"

namespace Cantera
{

struct FreezingData: public ReactionData
{
  FreezingData(): freezing_rate_(0.), temperature_(0.)
  {}

  virtual void update(double T) override;

  virtual void update(double T, double P) override;

  virtual bool update(const ThermoPhase& phase, const Kinetics& kin) override;

  using ReactionData::update;

protected:
  double freezing_rate_;
  double temperature_;
};

class FreezingRate final : public ReactionRate
{
public:
  FreezingRate();

  FreezingRate(const AnyMap& node, const UnitStack& rate_units={}) : 
    FreezingRate()
  {
      setParameters(node, rate_units);
  }

  unique_ptr<MultiRateBase> newMultiRate() const {
    return unique_ptr<MultiRateBase>(new MultiRate<FreezingRate, FreezingData>);
  }

  //! Identifier of reaction rate type
  const std::string type() const { return "frezing"; }
  
  void setParameters(const AnyMap& node, const UnitStack& rate_units);

  void getParameters(AnyMap& rateNode, const Units& rate_units) const;
  void getParameters(AnyMap& rateNode) const {
      return getParameters(rateNode, Units(0));
  }

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

  //! Set up Freezing object
  //void setRates(const std::multimap<double, ArrheniusRate>& rates);

  //! Check to make sure that the rate expression is finite over a range of
  //! temperatures at each interpolation pressure. This is potentially an
  //! issue when one of the Arrhenius expressions at a particular pressure
  //! has a negative pre-exponential factor.
  void validate(const std::string& equation, const Kinetics& kin);
};

}

#endif
