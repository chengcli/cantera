/**
 * @file Photolysis.h
 */

#ifndef CT_PHOTOLYSIS_H
#define CT_PHOTOLYSIS_H

int locate(double const *xx, double x, int n);

void interpn(double *val, double const *coor, double const *data,
             double const *axis, size_t const *len, int ndim, int nval);

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
    IABSORB = 0,
    IDISSOC = 1,
    IIONIZE = 2
  };

  PhotolysisData();

  bool update(const ThermoPhase& thermo, const Kinetics& kin) override;
  using ReactionData::update;

  //! \brief wavelength grid
  //!
  //! The wavelength grid is a vector of size nwave.
  //! Default units are nanometers.
  vector<double> wavelength;

  //! \brief actinic flux
  //!
  //! The actinic flux is a vector of size nwave.
  //! Default units are photons cm^-2 s^-1 nm^-1.
  vector<double> actinicFlux;
};

class PhotolysisBase : public ReactionRate {
 public:
  //! Default constructor
  PhotolysisBase() {}

  //! Constructor.
  /*!
   * @param temp Temperature grid
   * @param wavelength Wavelength grid
   * @param cross_section Cross-section data
   */
  PhotolysisBase(vector<double> const& temp, vector<double> const& wavelength,
                 vector<vector<double>> const& cross_section);

  //! Constructor based on AnyValue content
  PhotolysisBase(AnyValue const& rate, UnitSystem& units,
                 UnitStack& rate_units);

  explicit PhotolysisBase(AnyMap const& node, UnitStack const& rate_units={});

  void setParameters(AnyMap const& node, UnitStack const& rate_units) override;

  void getParameters(AnyMap& node) const override;

  void check(string const& equation) override;

  void validate(const string& equation, const Kinetics& kin) override;

 protected:
  //! number of temperature grid points
  size_t m_ntemp;

  //! number of wavelength grid points
  size_t m_nwave;

  //! temperature grid followed by wavelength grid
  vector<double> m_temp_wave_grid;

  //! \brief photolysis cross-section data
  //!
  //! The cross-section data is a two dimensional table of size (ntemp, nwave, 3).
  //! Each row contains the wavelength, photo-absorption cross-section
  //! photo-dissociation cross-section, and photo-ionization cross-section.
  //! Default units are nanometers, cm^2, cm^2, and cm^2, respectively.
  vector<double> m_crossSection;
};

//! Photolysis reaction rate type depends on temperature and the actinic flux
/*! 
 * A reaction rate coefficient of the following form.
 *
 * \f[
 *    k(T) = \int_{\lambda_1}^{\lambda_2} \sigma(\lambda) \phi(\lambda) d\lambda
 * \f]
 *
 * where \f$ \sigma(\lambda) \f$ is the cross-section and \f$ \phi(\lambda) \f$
 * is the actinic flux. \f$ \lambda_1 \f$ and \f$ \lambda_2 \f$ are the lower
 * and upper bounds of the wavelength grid.
 */
class PhotolysisRate : public PhotolysisBase {
 public:
  using PhotolysisBase::PhotolysisBase;  // inherit constructor

  unique_ptr<MultiRateBase> newMultiRate() const override {
    return make_unique<MultiRate<PhotolysisRate, PhotolysisData>>();
  }

  const string type() const override {
    return "Photolysis";
  }

  double evalFromStruct(PhotolysisData const& data) const {
    double wmin = m_temp_wave_grid[m_ntemp];
    double wmax = m_temp_wave_grid.back();

    int iwmin = locate(data.wavelength.data(), wmin, data.wavelength.size());
    int iwmax = locate(data.wavelength.data(), wmax, data.wavelength.size());

    double cross1, cross2;
    double coord[2] = {data.temperature, data.actinicFlux[iwmin]};
    size_t len[2] = {m_ntemp, m_nwave};

    interpn(&cross1, coord, m_crossSection.data(), m_temp_wave_grid.data(), len, 2, 1);

    double rate = 0.;
    for (int i = iwmin; i < iwmax; i++) {
      coord[1] = data.actinicFlux[i+1];
      interpn(&cross2, coord, m_crossSection.data(), m_temp_wave_grid.data(), len, 2, 1);
      rate += 0.5 * (cross1 * data.actinicFlux[i] + cross2 * data.actinicFlux[i+1]) 
                  * (data.wavelength[i+1] - data.wavelength[i]);
      cross1 = cross2;
    }

    return rate;
  }
};

}

#endif  // CT_PHOTOLYSIS_H
