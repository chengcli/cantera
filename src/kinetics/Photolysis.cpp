//! @file Photolysis.cpp

#include "cantera/kinetics/Photolysis.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/stringUtils.h"

namespace Cantera
{

bool PhotolysisData::update(const ThermoPhase& thermo, const Kinetics& kin)
{
  bool changed = false;
  double T = thermo.temperature();
  if (T != temperature) {
    update(T);
    changed = true;
  }

  return changed;
}

bool PhotolysisData::check() const
{
    // Check that the wavelength grid is valid
    if (wavelength.empty()) {
        throw CanteraError("PhotolysisData::update",
                           "Wavelength grid is empty.");
    }

    if (wavelength.size() < 2) {
        throw CanteraError("PhotolysisData::update",
                           "Wavelength grid must have at least two points.");
    }

    if (wavelength[0] <= 0.0) {
        throw CanteraError("PhotolysisData::update",
                           "Wavelength grid must be positive.");
    }

    for (size_t i = 1; i < wavelength.size(); i++) {
        if (wavelength[i] <= wavelength[i-1]) {
            throw CanteraError("PhotolysisData::update",
                               "Wavelength grid must be strictly increasing.");
        }
    }

    // Check that the actinic flux is valid
    if (actinicFlux.empty()) {
        throw CanteraError("PhotolysisData::update",
                           "Actinic flux is empty.");
    }

    if (actinicFlux.size() != wavelength.size()) {
        throw CanteraError("PhotolysisData::update",
                           "Actinic flux must have the same size as the wavelength grid.");
    }

    for (size_t i = 0; i < actinicFlux.size(); i++) {
        if (actinicFlux[i] < 0.0) {
            throw CanteraError("PhotolysisData::update",
                               "Actinic flux must be non-negative.");
        }
    }

    return true;
}

PhotolysisBase::PhotolysisBase(
    vector<double> const& temp, vector<double> const& wavelength,
    vector<std::string> const& branches,
    vector<double> const& cross_section):
  m_crossSection(cross_section)
{
  m_ntemp = temp.size();
  m_nwave = wavelength.size();

  m_temp_wave_grid.resize(m_ntemp + m_nwave);
  for (size_t i = 0; i < m_ntemp; i++) {
    m_temp_wave_grid[i] = temp[i];
  }

  for (size_t i = 0; i < m_nwave; i++) {
    m_temp_wave_grid[m_ntemp + i] = wavelength[i];
  }

  for (auto const& branch : branches) {
    m_branch.push_back(parseCompString(branch));
  }

  if (m_ntemp * m_nwave * branches.size() != m_crossSection.size()) {
    throw CanteraError("PhotolysisBase::PhotolysisBase",
                       "Cross-section data size does not match the temperature, "
                       "wavelength, and branch grid sizes.");
  }
}

PhotolysisBase::PhotolysisBase(AnyMap const& node, UnitStack const& rate_units)
{
  setParameters(node, rate_units);
}

void PhotolysisBase::setRateParameters(const AnyValue& rate, vector<string> const& branches)
{
  if (rate.hasKey("resolution")) {
    double resolution = rate["resolution"].asDouble();
    if (resolution <= 0.0) {
      throw CanteraError("PhotolysisBase::setRateParameters",
                         "Resolution must be positive.");
    }
  }

  if (rate.hasKey("scale")) {
    vector<double> scales;
    if (rate["scale"].is<double>()) {
      for (size_t i = 0; i < m_branch.size(); i++) {
        scales[i] = rate["scale"].asDouble();
      }
    } else {
      for (auto const& [branch, scale] : rate["scale"].as<AnyMap>()) {
        auto b = std::find(branches.begin(), branches.end(), branch);
        if (b == branches.end()) {
          throw CanteraError("PhotolysisBase::setRateParameters",
                             "Branch '{}' not found", branch);
        }

        scales[b - branches.begin()] = scale.asDouble();
      }
    }
  }
}

void PhotolysisBase::setParameters(AnyMap const& node, UnitStack const& rate_units)
{
  vector<string> branches;
  vector<double> temperature;
  vector<double> wavelength;
  vector<double> cross_section;

  ReactionRate::setParameters(node, rate_units);

  if (node.hasKey("branches")) {
    for (auto const& branch : node["branches"].asVector<AnyMap>()) {
      branches.push_back(branch["name"].asString());
      m_branch.push_back(parseCompString(branch["product"].asString()));
    } 
  } else {
    branches.push_back("all");
    vector<string> tokens;
    tokenizeString(node["equation"].asString(), tokens);
    m_branch.push_back(parseCompString(tokens[0]));
  }

  if (node.hasKey("cross-section")) {
    for (auto const& data: node["cross-section"].asVector<AnyValue>()) {
      auto format = data["format"].asString();
      std::string branch = data.hasKey("branch") ? data["branch"].asString() : "all";

      if (format == "YAML") {
        for (auto const& entry: data["data"].asVector<AnyValue>()) {
          auto values = entry.asVector<double>(2);
          wavelength.push_back(values[0]);
          cross_section.push_back(values[1]);
        }
      } else if (format == "VULCAN") {
        auto files = data["filenames"].asVector<string>();
        loadCrossSectionVulcan(files, branch);
      } else if (format == "KINETICS7") {
        auto files = data["filenames"].asVector<string>();
        loadCrossSectionKinetics7(files, branch);
      } else {
        throw CanteraError("PhotolysisBase::setParameters",
                           "Invalid cross-section format '{}'.", format);
      }
    }
  }

  m_ntemp = temperature.size();
  m_nwave = wavelength.size();
  m_temp_wave_grid.resize(m_ntemp + m_nwave);
  m_crossSection = cross_section;

  if (node.hasKey("rate-constant")) {
    setRateParameters(node["rate-constant"], branches);
  }
}

}
