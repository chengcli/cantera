/**
 * @file PhotolysisData.h
 */

#ifndef CT_PHOTOLYSIS_H
#define CT_PHOTOLYSIS_H

namespace Cantera
{

class ThermoPhase;
class Kinetics;

//! Data container holding shared data specific to photolysis reactions
/**
 * The data container `PhotolysisData` holds photolysis cross-section data
 * @ingroup reactionGroup
 */
struct PhotolysisData : public ReactionData {
  enum {
    IWAVE = 0,
    IABSORB = 1,
    IDISSOC = 2,
    IIONIZE = 3
  }

  PhotolysisData()

  bool update(const ThermoPhase& thermo, const Kinetics& kin) override;
  using ReactionData::update;

 protected:
  //! \brief photolysis cross-section data
  //!
  //! The cross-section data is a two dimensional table of size (nwave, 4).
  //! Each row contains the wavelength, photo-absorption cross-section
  //! photo-dissociation cross-section, and photo-ionization cross-section.
  vector<double> m_crossSection;
};

class PhotolysisBase : public ReactionRate {
 public:
  PhotolysisBase() {}
};

class PhotolysisRate : public PhotolysisBase {
 public:
  using PhotolysisBae::PhotolysisBase;

  unique_ptr<MultiRateBase> newMultiRate() const override {
    return make_unique<MultiRate<PhotolysisRate, PhotolysisData>>();
  }

  const string type() const override {
    return "Photolysis";
  }
};

}

#endif
